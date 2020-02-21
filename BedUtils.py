#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################
#     handle the bed format row data     #
#          2020.2.19                     #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.0"

import sys, os

class raiseErr(object):
    def __init__(self, message, exit=True):
        sys.stderr.write(message)
        if bool(exit):
            sys.exit()

class creatBed(object):
    def __init__(self, row):
        self.clear = True
        self.name, self.score, self.bcount, self.bsize, self.bstart = [None for i in range(5)]
        self.strand = '.'
        
        try:
            self.chr, self.start, self.end = row[0:3]
        except IndexError as e:
            self.clear = False
        
        if len(row) == 4:
            self.score = row[3]
        elif len(row) == 5:
            self.name, self.score = row[3:5]
        elif len(row) >= 6:
            self.name, self.score, self.strand = row[3:6]
        if len(row) >= 12:
            try:
                self.tstart = int(row[6])
                self.tend = int(row[7])
                self.bcount = int(row[9])
                self.bsize = [int(i) for i in row[10].strip(',').split(',') if i]
                self.bstart = [int(i) for i in row[11].strip(',').split(',') if i]
            except ValueError as e:
                self.clear = False

    # check bed
    def check(self):
        if self.clear:
            try:
                self.start = int(self.start)
            except (ValueError, TypeError) as e:
                self.clear = False
            
            try:
                self.end = int(self.end)
            except (ValueError, TypeError) as e:
                self.clear = False
            
            try:
                if self.score is not None:
                    self.score = float(self.score)
            except ValueError as e:
                self.clear = False
            
            try:
                if self.strand not in ['+', '-', '.']:
                    self.clear = False
                elif self.end <= self.start and (self.start < 0 or self.end <= 0):
                    self.clear = False
                elif self.bstart is not None and self.strand not in ['+', '-']:
                    self.clear = False
                else:
                    self.clear = True
            except TypeError as e:
                self.clear = False
        return self

    # decode bed12 to [ [exonblocks], [intronblocks] ] for ncRNA,
    # [ [exonblocks], [intronblocks], [ [5'utrs], [thicks], [3'utrs] ] ] for protein-coding
    # only return genomic interval in self.code
    # intronList will be empty if no intron
    # strand is consider
    def decode(self):
        ## used for test
        ## row = ['chr1','8423769','8424898','ENST00000464367','1000','-','8423770', '8424810','0','2','546,93,','0,1036,']
        def overlap(a, b):
            # 0-based
            distance = min(a[1], b[1]) - max(a[0], b[0])
            return max(0, distance)

        self.code = None
        if self.check().clear:
            blockList = list(map(lambda x,y:[x + self.start, x + self.start + y],
                self.bstart, self.bsize))
            intronList = list()
            if len(blockList) > 1:
                for i in range(len(blockList) - 1):
                    ## [exon.end, next.exon.start]
                    intronList.append([blockList[i][1], blockList[i+1][0]])
            if self.tstart == self.tend:
                if self.strand == '-':
                    blockList = blockList.reverse()
                    intronList = intronList.reverse()
                decodeList = [blockList, intronList]
            else:
                ## for protein-coding transcript
                decodeList = [blockList, intronList, [[], [], []]]
                thickStartLocus = [self.tstart, self.tstart + 1]
                thickEndLocus = [self.tend - 1, self.tend]
                ## return thick-start and thick-end exon index
                thickStartBoolIndex = list(map(lambda x:overlap(thickStartLocus, x), blockList)).index(1)
                thickEndBoolIndex = list(map(lambda x:overlap(thickEndLocus, x), blockList)).index(1)
                for i in range(len(blockList)):
                    blockStart = blockList[i][0]
                    blockEnd = blockList[i][1]
                    if i < thickStartBoolIndex:
                        decodeList[-1][0].append([blockStart, blockEnd])
                    elif i == thickStartBoolIndex:
                        if self.tstart > blockStart:
                            decodeList[-1][0].append([blockStart, self.tstart])
                        if i == thickEndBoolIndex:
                            decodeList[-1][1].append([self.tstart, self.tend])
                            if self.tend < blockEnd:
                                decodeList[-1][2].append([self.tend, blockEnd])
                        else:
                            decodeList[-1][1].append([self.tstart, blockEnd])
                    elif i > thickStartBoolIndex and i < thickEndBoolIndex:
                        decodeList[-1][1].append([blockStart, blockEnd])
                    elif i == thickEndBoolIndex:
                        decodeList[-1][1].append([blockStart, self.tend])
                        if self.tend < blockEnd:
                            decodeList[-1][2].append([self.tend, blockEnd])
                    else:
                        decodeList[-1][2].append([blockStart, blockEnd])
            if self.strand == '-':
                for i in range(2):
                    decodeList[i].reverse()
                if self.tstart != self.tend:
                    for i in range(3):
                        decodeList[2][i].reverse()
                    decodeList[2][0], decodeList[2][-1] = decodeList[2][-1], decodeList[2][0]
            self.code = decodeList
        return self

# bed operations on 2 bed6 format row
class bedops(object):
    # a:bed locus-A, b:bed locus-B
    # s:strand, d:distance
    def __init__(self, a, b, s=False):
        self.a = creatBed(a)
        self.b = creatBed(b)
        self.strand = s
        self.clear = True
    # check intervals
    def check(self):
        if self.a.check().clear and self.b.check().clear:
            if self.a.chr != self.b.chr:
                self.clear = False
            if self.strand:
                if self.a.strand != self.b.strand:
                    self.clear = False
        else:
            self.clear = False
        return self
    # calculate distance between intervals
    def discompute(self, tss=False, center=False):
        # return distance ralative to a
        # if tss (transcription start site) is set True, take locus-B as genomic locus, a as RNA-type locus
        # tss will ignore self.strand
        # center only work with tss
        self.distance = None
        if tss:
            self.strand = False
        if self.check().clear:
            overlapLength = self.intersect().overlap
            if overlapLength > 0:
                distance = 0
            else:
                if tss:
                    if center:
                        peak = int((self.b.end + self.b.start) / 2)
                        distance = (peak - self.a.end) if self.a.strand == '-' else (peak - self.a.start)
                    else:
                        if self.a.strand == '+' or self.a.strand == '.':
                            distance = min(abs(self.b.start - self.a.start), abs(self.b.end - self.a.start))
                            if self.b.end <= self.a.start: distance = -distance
                        else:
                            distance = min(abs(self.b.start - self.a.end + 1), abs(self.b.end - self.a.end))
                            if self.a.end <= self.b.start: distance = -distance
                else:
                    distance = min(abs(self.b.start - self.a.end + 1), abs(self.b.end - self.a.start))
                    if self.b.end <= self.a.start: distance = -distance
                    if self.strand and self.a.strand == '-': distance = -distance
            self.distance = distance
        return self
    # calcualte intersection length of intervals
    def intersect(self):
        # 0-based, return 0 if no overlap
        self.overlap = None
        if self.check().clear:
            length = min(self.a.end, self.b.end) - max(self.a.start, self.b.start)
            self.overlap = max(0, length)
            self.fracA = length / (self.a.end - self.a.start)
            self.fracB = length / (self.b.end - self.b.start)
            self.ilocus = None
            if self.overlap > 0:
                self.ilocus = [max(self.a.start, self.b.start), min(self.a.end, self.b.end)]
        return self
    # merge intervals
    def merge(self, d=0):
        self.locus = None
        if self.check().clear:
            overlapLength = self.intersect().overlap
            distance = abs(self.discompute().distance)
            if overlapLength > 0 or distance <= d:
                chrom = self.a.chr
                name = ':'.join([self.a.name, self.b.name])
                score = self.a.score + self.b.score
                if self.strand:
                    strand = self.a.strand
                else:
                    strand = '.'
                locus = [chrom, min(self.a.start, self.b.start), max(self.a.end, self.b.end), name, score, strand]
            else:
                locus = None
            self.locus = locus
        return self
    # calculate the information of how a include b (required intersections)
    def include(self):
        pass
