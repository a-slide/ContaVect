#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from multiprocessing import Value, Process, Queue, cpu_count
from time import time
import gzip
from os import path

# Third party package import
from Bio import SeqIO

# Local Package import
from pyDNA.Utilities import file_basename, count_seq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqFilter(object):
    """
    @class  FastqFilter
    @brief Main class of the package
    Require the third party package Biopython
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, R1, R2,
        quality_filter=None,
        adapter_trimmer=None,
        outdir="./fastq/",
        input_qual="fastq-sanger",
        numprocs=None,
        compress_output=True):
        """
        Instanciate the object by storing call parameters and init shared memory counters for
        interprocess communication. A reader process iterate over the input paired fastq files
        and add coupled R1 and R2 sequences as Biopython seqRecord to a first shared queue.
        Then according to the initial parametring, a multiprocessing filter pull out seqRecord
        couples from the queue and apply a quality filtering and/or adapter trimming. Couples
        passing throught the filters are added to a second shared queue. Finally, couples in the
        second queue are written in an output fastq file
        @param R1 Path to the forward read fastq file (can be gzipped)
        @param R2 Path to the reverse read fastq file (can be gzipped)
        @param quality_filter A QualityFilter object, if a quality filtering is required.
        @param adapter_trimmer An AdapterTrimmer object, if a adapter trimming is required.
        @param outdir Directory where to write the filtered fastq sequences.
        @param input_qual Quality scale of the fastq (fastq-sanger for illumina 1.8+)
        @param numprocs Number of parrallel processes for the filtering steps. If not provide
        the maximum number of thread available will be automatically used.
        @param compress_output If True the output fastq will be written directly in a gzipped file.
        False will generate an uncompressed a much bigger file but will be around
        """
        # Start a timer
        start_time = time()

        # Create object variables
        self.numprocs = numprocs if numprocs else cpu_count()
        self.qual = quality_filter
        self.adapt = adapter_trimmer
        self.input_qual = input_qual
        self.R1_in = R1
        self.R2_in = R2
        self.outdir = outdir
        self.compress_output = compress_output
        if compress_output:
            self.R1_out = path.join(self.outdir, file_basename(self.R1_in)+"_1_filtered.fastq.gz")
            self.R2_out = path.join(self.outdir, file_basename(self.R2_in)+"_2_filtered.fastq.gz")
        else:
            self.R1_out = path.join(self.outdir, file_basename(self.R1_in)+"_1_filtered.fastq")
            self.R2_out = path.join(self.outdir, file_basename(self.R2_in)+"_2_filtered.fastq")

        # Init shared memory counters
        self.total = Value('i', 0)
        self.pass_qual = Value('i', 0)
        self.pass_trim = Value('i', 0)
        self.total_pass = Value('i', 0)
        if self.qual:
            self.min_qual_found = Value('i', 100)
            self.max_qual_found = Value('i', 0)
            self.weighted_mean = Value('d', 0.0)
        if self.adapt:
            self.seq_untrimmed = Value('i', 0)
            self.seq_trimmed = Value('i', 0)
            self.base_trimmed = Value('i', 0)
            self.len_pass = Value('i', 0)
            self.len_fail = Value('i', 0)

        # Count lines in fastq file to prepare a counter of progression
        print ("Count the number of fastq sequences")
        self.nseq = count_seq(R1, "fastq")
        print(("fastq files contain {} sequences to align".format(self.nseq)))
        self.nseq_list = [int(self.nseq*i/100.0) for i in range(5,101,5)] # 5 percent steps

        # Init queues for input file reading and output file writing (limited to 10000 objects)
        self.inq = Queue(maxsize=10000)
        self.outq = Queue(maxsize=10000)

        # Init processes for file reading, distributed filtering and file writing
        self.pin = Process(target=self.reader, args=())
        self.ps = [Process(target=self.filter, args=()) for i in range(self.numprocs)]
        self.pout = Process(target=self.writer, args=())

        # Start processes
        self.pin.start()
        self.pout.start()
        for p in self.ps:
            p.start()

        # Blocks until the process is finished
        self.pin.join()
        print ("\tReading done")
        for i in range(len(self.ps)):
            self.ps[i].join()
        print ("\tFiltering done")
        self.pout.join()
        print ("\tWriting done\n")

        # Stop timer and store the value
        self.exec_time = round(time()-start_time, 3)

    def __repr__(self):
        msg = "FASTQ FILTER Parallel Processing\n"
        msg += "\tExecution time : {} s\n".format(self.exec_time)
        msg += "\tInput fastq files\n\t\t{}\n\t\t{}\n".format (self.R1_in, self.R2_in)
        msg += "\tOutput fastq files\n\t\t{}\n\t\t{}\n".format (self.R1_out, self.R2_out)
        msg += "\tInput quality score : {}\n".format (self.input_qual)
        msg += "\tNumber of parallel processes : {}\n".format (self.numprocs)
        msg += "\tTotal pair processed : {}\n".format(self.total.value)
        msg += "\tTotal pair passed : {}\n".format(self.total_pass.value)
        if self.qual:
            msg += "QUALITY FILTER\n"
            msg += "\tPair pass quality filter : {}\n".format(self.pass_qual.value)
            msg += "\tMean quality value : {}\n".format(self.weighted_mean.value/self.total.value/2)
            msg += "\tMin quality value : {}\n".format(self.min_qual_found.value)
            msg += "\tMax quality value : {}\n".format(self.max_qual_found.value)
        if self.adapt:
            msg += "ADAPTER TRIMMER\n"
            msg += "\tPair pass adapter Trimming : {}\n".format(self.pass_trim.value)
            msg += "\tSequences untrimmed : {}\n".format(self.seq_untrimmed.value)
            msg += "\tSequences trimmed : {}\n".format(self.seq_trimmed.value)
            msg += "\tDNA base trimmed : {}\n".format(self.base_trimmed.value)
            msg += "\tFail len filtering: {}\n".format(self.len_fail.value)
            msg += "\tPass len filtering : {}\n".format(self.len_pass.value)
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def get(self, key):
        return self.__dict__[key]

    def getCTypeVal(self, key):
        return self.__dict__[key].value
    
    def getTrimmed (self):
        return (self.R1_out, self.R2_out)

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def reader(self):
        """
        Initialize SeqIO.parse generators to iterate over paired fastq files. Data ara sent over
        inqueue for the workers to do their thing and a n = numprocs STOP pills are added at the
        end of the queue for each worker.
        """
        try:
            # Open input fastq streams for reading
            if self.R1_in[-2:].lower() == "gz":
                in_R1 = gzip.open(self.R1_in, "r")
            else:
                in_R1 = open(self.R1_in, "r")

            if self.R2_in[-2:].lower() == "gz":
                in_R2 = gzip.open(self.R2_in, "r")
            else:
                in_R2 = open(self.R2_in, "r")

        except (IOError, TypeError, ValueError) as E:
            print(E)
            exit

        # Init generators to iterate over files
        genR1 = SeqIO.parse(in_R1, self.input_qual)
        genR2 = SeqIO.parse(in_R2, self.input_qual)

        i = 0
        while True:
            # Parse sequences in generators until one of then is empty
            seqR1 = next(genR1, None)
            seqR2 = next(genR2, None)
            if not seqR1 or not seqR2:
                break
            # Add a tuple position, seqR1 and seqR2 to the end of the queue
            self.inq.put( (seqR1, seqR2) )

            i+=1
            if i in self.nseq_list:
                print(("\t{} sequences: {}%".format(i, int(i*100.0/self.nseq))))

        # Close files
        in_R1.close()
        in_R2.close()

        # Add a STOP pill to the queue
        for i in range(self.numprocs):
            self.inq.put("STOP")

    def filter(self):
        """
        Parallelized filter that take as input a sequence couple in inqueue until a STOP pill is
        found. Sequences go through a QualityFilter and a AdapterTrimmer object and ifthe couple
        is able to pass filters then it is put at the end of outqueue. at the ebd of the process
        a STOP pill is added to the outqueue.
        """
        # Consume inq and produce answers on outq
        for seqR1, seqR2 in iter(self.inq.get, "STOP"):

            with self.total.get_lock():
                self.total.value+=1

            # Quality filtering
            if self.qual:
                seqR1 = self.qual.filter(seqR1)
                seqR2 = self.qual.filter(seqR2)
                if not seqR1 or not seqR2:
                    continue

            with self.pass_qual.get_lock():
                self.pass_qual.value+=1

            # Adapter trimming and size filtering
            if self.adapt:
                seqR1 = self.adapt.trimmer(seqR1)
                seqR2 = self.adapt.trimmer(seqR2)
                if not seqR1 or not seqR2:
                    continue

            with self.pass_trim.get_lock():
                self.pass_trim.value+=1

            # If both filters passed = add to the output queue
            self.outq.put( (seqR1, seqR2) )

        # Add a STOP pill to the queue
        self.outq.put("STOP")

        # Fill shared memory counters from process specific object instances.
        if len(self.qual.mean_qual) > 0:
            with self.weighted_mean.get_lock():
                self.weighted_mean.value += (self.qual.get_mean_qual()*self.qual.get('total'))
            if self.qual.get_min_qual() < self.min_qual_found.value:
                self.min_qual_found.value = int(self.qual.get_min_qual())
            if self.qual.get_max_qual() > self.max_qual_found.value:
                self.max_qual_found.value = int(self.qual.get_max_qual())

        if self.adapt:
            with self.seq_untrimmed.get_lock():
                self.seq_untrimmed.value += self.adapt.get('seq_untrimmed')
            with self.seq_trimmed.get_lock():
                self.seq_trimmed.value += self.adapt.get('seq_trimmed')
            with self.base_trimmed.get_lock():
                self.base_trimmed.value += self.adapt.get('base_trimmed')
            with self.len_pass.get_lock():
                self.len_pass.value += self.adapt.get('len_pass')
            with self.len_fail.get_lock():
                self.len_fail.value += self.adapt.get('len_fail')

    def writer(self):
        """
        Write sequence couples from outqueue in a pair of compressed fastq.gz files. Sequences will
        remains paired (ie at the same index in the 2 files) but they may not be in the same order
        than in the input fastq files. The process will continue until n = numprocs STOP pills were
        found in the outqueue (ie. the queue is empty)
        """
        # Open output fastq streams for writing
        if self.compress_output:
            out_R1 = gzip.open(self.R1_out, "wb")
            out_R2 = gzip.open(self.R2_out, "wb")
        else:
            out_R1 = open(self.R1_out, "wb")
            out_R2 = open(self.R2_out, "wb")

        # Keep running until all numprocs STOP pills has been passed
        for works in range(self.numprocs):
            # Will exit the loop as soon as a Stop pill will be found
            for seqR1, seqR2 in iter(self.outq.get, "STOP"):
                out_R1.write(seqR1.format("fastq-sanger").encode())
                out_R2.write(seqR2.format("fastq-sanger").encode())
                with self.total_pass.get_lock():
                    self.total_pass.value+=1

        out_R1.close()
        out_R2.close()

# Required by multiprocessing
if __name__ == '__main__':
    pass
