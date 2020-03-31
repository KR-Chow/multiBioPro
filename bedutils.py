# -*- coding: utf-8 -*-
##########################################
#     handle the bed format row data     #
#          2020.2.19                     #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.0"

import random

class creatbed(object):
    def __init__(self, row):
        self.clear = True
        self.number = len(row)
        self.name = str(random.randrange(100000))
        self.score = 0
        self.bcount, self.bsize, self.bstart = [None for i in range(3)]
        self.strand = '.'
        try:
            self.chr, self.start, self.end = row[0:3]
        except IndexError as e:
            self.clear = False
        if self.number == 4:
            self.score = row[3]
        elif self.number == 5:
            self.name, self.score = row[3:5]
        elif self.number >= 6:
            self.name, self.score, self.strand = row[3:6]
        if self.number >= 12:
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
                self.end = int(self.end)
                self.length = self.end - self.start
                if self.score is not None:
                    self.score = float(self.score)
            except (ValueError, TypeError) as e:
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
        if self.clear:
            return self
        else:
            raise SystemExit("Error when passing row! Please pass bed-like row to creatbed!")
        return self

    # strand is taken into consideration
    # if "-" strand, all coordinates were re-order by reversing
    # row = ['chr1','8423769','8424898','ENST00000464367','1000','-','8423770', '8424810','0','2','546,93,','0,1036,']
    # self.exon = [[8424805, 8424898], [8423769, 8424315]]
    # self.exon = [[8424315, 8424805]]
    # self.exon = [[8424805, 8424810], [8423770, 8424315]]
    # self.exon = [[8424810, 8424898]]
    # self.exon = [[8423769, 8423770]]
    def decode(self):
        ## return overlap length
        def overlap(a, b):
            # 0-based
            distance = min(a[1], b[1]) - max(a[0], b[0])
            return max(0, distance)
        ## main code
        self.exon = list()
        self.intron = list()
        self.cds = list()
        self.utr5 = list()
        self.utr3 = list()
        if self.check().clear:
            self.exon = list(map(lambda x,y:[x + self.start, x + self.start + y],
                self.bstart, self.bsize))
            self.intron = list()
            if len(self.exon) > 1:
                for i in range(len(self.exon) - 1):
                    ## [exon.end, next.exon.start]
                    self.intron.append([self.exon[i][1], self.exon[i+1][0]])
            ## thick start and thick end
            if self.tstart == self.tend:
                ## for non-coding transcript
                if self.strand == '-':
                    self.exon = self.exon.reverse()
                    self.intron = self.intron.reverse()
            else:
                ## for protein-coding transcript
                tstartLocus = [self.tstart, self.tstart + 1]
                tendLocus = [self.tend - 1, self.tend]
                ## return exon index of thick-start and thick-end
                tstartExon = list(map(lambda x:overlap(tstartLocus, x), self.exon)).index(1)
                tendExon = list(map(lambda x:overlap(tendLocus, x), self.exon)).index(1)
                for i in range(len(self.exon)):
                    blockStart = self.exon[i][0]
                    blockEnd = self.exon[i][1]
                    if i < tstartExon:
                        self.utr5.append([blockStart, blockEnd])
                    elif i == tstartExon:
                        if self.tstart > blockStart:
                            self.utr5.append([blockStart, self.tstart])
                        if i == tendExon:
                            self.cds.append([self.tstart, self.tend])
                            if self.tend < blockEnd:
                                self.utr3.append([self.tend, blockEnd])
                        else:
                            self.cds.append([self.tstart, blockEnd])
                    elif i > tstartExon and i < tendExon:
                        self.cds.append([blockStart, blockEnd])
                    elif i == tendExon:
                        self.cds.append([blockStart, self.tend])
                        if self.tend < blockEnd:
                            self.utr3.append([self.tend, blockEnd])
                    else:
                        self.utr3.append([blockStart, blockEnd])
            if self.strand == '-':
                self.utr5, self.utr3 = self.utr3, self.utr5
                tempList = [self.exon, self.intron, self.cds, self.utr5, self.utr3]
                for coord in tempList:
                    coord.reverse()
        return self

# bed operations on 2 bed6 format row
class bedops(object):
    # a:bed locus-A, b:bed locus-B
    # s:strand, d:distance
    def __init__(self, a, b, s=False):
        self.a = creatbed(a)
        self.b = creatbed(b)
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
        if self.clear:
            return self
        else:
            emessage = "Error when passing rows!\nHere are some reasons:\n1.Not bed format;\n2.Not on the same chromosome;\n"
            emessage += "3.Not the same strand when s=True)"
            raise SystemExit(emessage)

    # calculate distance between intervals
    def discompute(self, tss=False, center=False):
        # tss=False, return distance ralative to genome (b to a), ignore strand
        # tss=True (transcription start site), take locus-B as genomic locus, a as RNA-type locus
        # tss will ignore self.strand
        # center only work with tss
        tssFlag, centerFlag = tss, center
        self.distance = None
        self.strand = False
        if self.check().clear:
            overlapLength = self.intersect().ilength
            if overlapLength > 0:
                distance = 0
            else:
                if tssFlag:
                    if centerFlag:
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
            self.distance = distance
        return self
    # calcualte intersection length of intervals
    def intersect(self):
        # 0-based, return 0 if no overlap
        # sensitive to strand
        self.ilength, self.ilocus, self.ifracA, self.ifracB = [None for i in range(4)]
        self.ibool = False
        if self.check().clear:
            length = min(self.a.end, self.b.end) - max(self.a.start, self.b.start)
            self.ilength = max(0, length)
            if self.ilength > 0:
                self.ibool = True
                self.ifracA = length / (self.a.end - self.a.start)
                self.ifracB = length / (self.b.end - self.b.start)
                name = ':'.join([self.a.name, self.b.name])
                score = max(self.a.score, self.b.score)
                if self.strand:
                    strand = self.a.strand
                else:
                    strand = '.'
                self.ilocus = [self.a.chr, max(self.a.start, self.b.start), min(self.a.end, self.b.end), name, score, strand]
        return self
    # merge intervals
    def merge(self, d=0):
        # if s=False, return strand with '.' when a and b has different orientations
        setDistance = d
        self.mlocus = None
        if self.check().clear:
            overlapLength = self.intersect().ilength
            distance = abs(self.discompute().distance)
            if overlapLength > 0 or distance <= setDistance:
                name = ':'.join([self.a.name, self.b.name])
                score = self.a.score + self.b.score
                if self.strand:
                    strand = self.a.strand
                else:
                    strand = '.'
                mlocus = [self.a.chr, min(self.a.start, self.b.start), max(self.a.end, self.b.end), name, score, strand]
            else:
                mlocus = None
            self.mlocus = mlocus
        return self
    # calculate the information of how a include b (required intersections)
    def include(self):
        # calculate how a include b
        # self.strand not work
        #cloverh: left overhang of locusB, croverh: right overhang of locusB
        # ctype:0->complete, 1->right, -1->left, 2->overlay
        self.strand = False
        self.cloverh = None
        self.croverh = None
        self.ctype = None
        if self.check().clear:
            intersectObj = self.intersect()
            overlapLength = intersectObj.ilength
            locusbLength = self.b.end - self.b.start
            if overlapLength > 0:
                cloverh = self.b.start - self.a.start
                croverh = self.b.end - self.a.end
                if cloverh >= 0:
                    if croverh <= 0:
                        ctype = 0
                    else:
                        ctype = 1
                else:
                    if croverh < 0:
                        ctype = -1
                    else:
                        ctype = 2
                self.cloverh = cloverh
                self.croverh = croverh
                self.ctype = ctype
        return self
