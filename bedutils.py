# -*- coding: utf-8 -*-
##########################################
#     handle the bed format row data     #
#          2021.6.5                      #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.0"

import random

class initbed(object):
    def __init__(self, n=6):
        self.chr = None
        self.start = None
        self.end = None
        self.name = str(random.randrange(100))
        self.score = 255
        self.strand = '.'
        if n > 6:
            self.tstart = None
            self.tend = None
            self.rgb = '255,0,0'
            self.bcount = None
            self.bsize = None
            self.bstart = None

class buildbed(object):
    def __init__(self, row):
        self.clear = True
        self.list = list(map(str, row))
        self.colnum = len(row)
        self.name = str(random.randrange(100))
        self.score = 255
        self.strand = '.'
        self.bcount = None
        self.bsize = None
        self.bstart = None
        try:
            self.chr, self.start, self.end = row[0:3]
        except IndexError as e:
            self.clear = False
        if self.colnum == 4:
            self.score = row[3]
        elif self.colnum == 5:
            self.name, self.score = row[3:5]
        elif self.colnum >= 6:
            self.name, self.score, self.strand = row[3:6]
        if self.colnum >= 12:
            try:
                self.tstart = int(row[6])
                self.tend = int(row[7])
                self.rgb = row[8]
                self.bcount = int(row[9])
                self.bsize = [int(i) for i in row[10].strip(',').split(',') if i]
                self.bstart = [int(i) for i in row[11].strip(',').split(',') if i]
            except ValueError as e:
                self.clear = False
        try:
            self.start = int(self.start)
            self.end = int(self.end)
            self.score =float(self.score)
        except ValueError as e:
            raise SystemError("The chromStart or chromEnd should be integer, or score should be number!")
       # check bed
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
            if self.colnum >= 12:
                if self.bcount != len(self.bsize) or self.bcount != len(self.bstart):
                    raise SystemError("Input is not a standard bed12 row!")
                if self.tstart < self.start or self.tend > self.end or self.tstart == self.bsize[0]:
                    raise SystemError("The thick start-end should not exceed the bed region!")
                self.exonlength = sum(self.bsize)
        if self.clear is False:
            raise SystemError("Error when passing row! Please pass bed-like row to buildbed!")
    # return coordinates are in bed format
    # row = ['chr1','8423769','8424898','ENST00000464367','1000','-','8423770', '8424810','0','2','546,93,','0,1036,']
    # self.exon = [[8423769, 8424315], [8424805, 8424898]]
    # self.intron = [[8424315, 8424805]]
    # self.cds = [[8423770, 8424315], [8424805, 8424810]]
    # self.utr5 = [[8424810, 8424898]]
    # self.utr3 = [[8423769, 8423770]]
    def decode(self):
        ## return overlap length
        def overlap(a, b):
            distance = min(a[1], b[1]) - max(a[0], b[0])
            if distance > 0:
                return True
            else:
                return False
        ## main code
        self.exon = []
        self.intron = []
        self.cds = []
        self.utr5 = []
        self.utr3 = []
        self.exon = list(map(lambda x,y:[x + self.start, x + self.start + y], self.bstart, self.bsize))
        if self.exon[0][0] != self.start or self.exon[-1][-1] != self.end:
            raise SystemError("Input is not a standard bed12 row!")
        ## check if there are overlaps between exons
        blockDict = {}
        for exon in self.exon:
            for i in range(exon[0], exon[1] + 1):
                if i in blockDict:
                    raise SystemError("There are overlapped exons in input bed12 row!")
                else:
                    blockDict[i] = 1
        self.intron = []
        if len(self.exon) > 1:
            for i in range(len(self.exon) - 1):
                ## [exon.end, next.exon.start]
                self.intron.append([self.exon[i][1], self.exon[i+1][0]])
        ## thick start and thick end
        if self.tstart != self.tend:
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
        return self

# bed operations on 2 bed6 format row
class bed6ops(object):
    # a:bed locus-A, b:bed locus-B
    # s:strand, d:distance
    # retrun self.a self.b, self.i, self.m
    def __init__(self, a):
        self.a = buildbed(a)
        self.list = self.a.list
        self.clear = self.a.clear
        if self.clear is False:
            raise SystemError("Input parameters could not pass the requirment!")
    ## check bed
    def __check(self):
        if self.strand is not True and self.strand is not False:
            raise SystemError("Input parameters 's' should be bool type.")
        else:
            self.clear = True
        ## check bed12 for self.a and self.b
        if self.b.clear:
            if self.a.chr != self.b.chr:
                self.clear = False
            elif self.strand:
                if self.a.strand != self.b.strand:
                    self.clear = False
        else:
            self.clear = False
        ## throw erros if inputs are not in bed12 format
        if self.clear is False:
            raise SystemError("Input parameters could not pass the requirment!")
        return self.clear
    # calculated score
    def __getScore(self, scoreA, scoreB, method='sum'):
        socreList = [scoreA, scoreB]
        if method == 'sum':
            score = sum(socreList)
        elif method == 'min':
            score = min(socreList)
        elif method == 'max':
            score = max(socreList)
        elif method == 'average':
            score = sum(socreList) / 2
        return score
    # return True if a overlap with b
    def __overlap(self):
        distance = min(self.a.end, self.b.end) - max(self.a.start, self.b.start)
        if distance > 0:
            return True
        else:
            return False
    # covert bed6 to bed12
    def tobed12(self, rgb='255,0,0'):
        tstart = self.a.start
        tend = self.a.end
        bcount = 1
        bsize = str(self.a.length) + ','
        bstart = '0,'
        bedrow = [self.a.chr, self.a.start, self.a.end, self.a.name, self.a.score, self.a.strand]
        bedrow += [tstart, tend, rgb, bcount, bsize, bstart]
        return buildbed(bedrow)
    # calculate distance between intervals
    def discompute(self, b, tss=False, center=False):
        # tss=False, return distance ralative to genome (b to a), ignore strand
        # tss=True (transcription start site), take locus-B as genomic locus, a as RNA-type locus
        # tss will ignore self.strand
        # center only work with tss
        # if tss:False, center:False, [1, 2] and [2, 4] returns 1
        # if tss:False, center:False, [1, 3] and [2, 4] returns 0
        # if tss:False, center:False, [1, 3] and [0, 1] returns -1
        self.b = buildbed(b)
        self.clear = self.__check()
        tssFlag, centerFlag = tss, center
        self.distance = None
        self.strand = False
        overlap = self.intersect(b)
        if bool(overlap) is True:
            overlapLength = overlap.ilength
        else:
            overlapLength = 0
        if overlapLength > 0:
            distance = 0
        else:
            if tssFlag:
                if centerFlag:
                    peak = int((self.b.end + self.b.start) / 2)
                    distance = (peak - self.a.end) if self.a.strand == '-' else (peak - self.a.start - 1)
                else:
                    if self.a.strand == '+' or self.a.strand == '.':
                        distance = min(abs(self.b.start - self.a.end + 1), abs(self.b.end - self.a.start - 1))
                        if self.b.end <= self.a.start:
                            distance = -distance
                    else:
                        distance = min(abs(self.b.start - self.a.end + 1), abs(self.b.end - self.a.end))
                        if self.a.end <= self.b.start:
                            distance = -distance
            else:
                distance = min(abs(self.b.start - self.a.end + 1), abs(self.b.end - self.a.start - 1))
                if self.b.end <= self.a.start:
                    distance = -distance
        self.distance = distance
        return self
    # calcualte intersection length of intervals
    def intersect(self, b, s=False, score='average'):
        ## 0-based, return False if no intersect
        ## sensitive to strand
        ## if s=False, return strand with '.' when a and b has different orientations
        self.b = buildbed(b)
        self.strand = s
        self.clear = self.__check()
        name = '|'.join([self.a.name, self.b.name])
        length = min(self.a.end, self.b.end) - max(self.a.start, self.b.start)
        length = max(0, length)
        if length > 0:
            fracA = length / (self.a.end - self.a.start)
            fracB = length / (self.b.end - self.b.start)
            newScore = self.__getScore(self.a.score, self.b.score, method=score)
            if self.strand:
                strand = self.a.strand
            else:
                strand = '.'
            chrom = self.a.chr
            start = max(self.a.start, self.b.start)
            end = min(self.a.end, self.b.end)
            row = [chrom, start, end, name, newScore, strand]
            self = bed6ops(row)
            ## indicate the function
            self.fun = 'intersect'
            self.ilength = length
            self.fracA = fracA
            self.fracB = fracB
            return self
        else:
            return False
    # merge intervals
    def merge(self, b, s=False, d=0, score='sum'):
        ## 0-based, return False if no merge
        ## always merge the overlapped blocks
        ## d=0, eg. [1,2], [2,3]
        ## d=-1, for overlap with overlapped blocks only
        ## if s=False, return strand with '.' when a and b has different orientations
        self.b = buildbed(b)
        self.strand = s
        self.clear = self.__check()
        setDistance = d
        overlap = self.__overlap()
        distance = abs(self.discompute(b).distance) - 1
        if overlap is True or distance <= setDistance:
            name = '|'.join([self.a.name, self.b.name])
            newScore = self.__getScore(self.a.score, self.b.score, method=score)
            if self.strand:
                strand = self.a.strand
            else:
                strand = '.'
            chrom = self.a.chr
            start = min(self.a.start, self.b.start)
            end = max(self.a.end, self.b.end)
            row = [chrom, start, end, name, newScore, strand]
            self = bed6ops(row)
            ## indicate the function
            self.fun = 'merge'
            return self
        else:
            return False
    # calculate the information of how a include b (required intersections)
    def include(self, b, s=False):
        # calculate how a include b
        # self.strand not work
        #cloverh: left overhang of locusB, croverh: right overhang of locusB
        # ctype:0->complete, 1->right, -1->left, 2->overlay
        self.b = buildbed(b)
        self.strand = s
        self.clear = self.__check()
        self.strand = False
        self.cloverh = None
        self.croverh = None
        self.ctype = None
        overlapLength = min(self.a.end, self.b.end) - max(self.a.start, self.b.start)
        overlapLength = max(0, overlapLength)
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

# bed operations on 2 bed12 format row from bedtools intersect -split -wa -wb -s
class bed12ops(object):
    # a:bed locus-A, b:bed locus-B
    # s:strand, d:distance
    # retrun self.a self.b, self.i, self.m
    def __init__(self, a):
        self.a = buildbed(a)
        self.list = self.a.list
        ## check bed12 for self.a
        self.clear = True
        if self.a.clear:
            if self.a.colnum != 12:
                self.clear = False
        else:
            self.clear = False
        ## throw erros if inputs are not in bed12 format
        if self.clear is False:
            raise SystemError("Input parameters could not pass the requirment!")
        ## get structures for self.a and self.b: exon, intron, cds, utr5, utr3
        self.a = self.a.decode()
    # check bed
    def __check(self):
        if self.strand is not True and self.strand is not False:
            raise SystemError("Input parameters 's' should be bool type.")
        else:
            self.clear = True
        ## check bed12 for self.a and self.b
        if self.b.clear:
            if self.a.chr != self.b.chr:
                self.clear = False
            elif self.strand:
                if self.a.strand != self.b.strand:
                    self.clear = False
            elif self.b.colnum != 12:
                self.clear = False
        else:
            self.clear = False
        ## throw erros if inputs are not in bed12 format
        if self.clear is False:
            raise SystemError("Input parameters could not pass the requirment!")
        return self.clear
    # get merged or overlap exons
    def __squeezeBlock(self, block1List, block2List, btype):
        ## store all block coordinates in a dictionary, ther value repsents the overlap counts
        if self.__tx is False:
            ## if not in tx mode, means we can freely suqeeze block
            blockDict = {}
            for block in block1List + block2List:
                for j in range(block[0], block[1] + 1):
                    if j not in blockDict:
                        blockDict[j] = 1
                    else:
                        blockDict[j] += 1
            ## construct the merged or overlapped exons
            if self.__fun == 'merge':
                positionList = sorted(blockDict.keys())
            elif self.__fun == 'intersect':
                positionList = sorted(filter(lambda x:blockDict[x] > 1, blockDict.keys()))
            if self.__overlap is not True and self.__overlap is not False:
                raise SystemError("Input parameters 'overlap' should be bool type.")
            if self.__tx is not True and self.__tx is not False:
                raise SystemError("Input parameters 'tx' should be bool type.")
            if self.__part is not True and self.__part is not False:
                raise SystemError("Input parameters 'part' should be bool type.")
            ## raise errer when no overlaps found restricted by 'overlap'
            if self.__overlap and btype == 'exon':
                if self.__fun != 'intersect':
                    positionList = list(filter(lambda x:blockDict[x] > 1, blockDict.keys()))
                if len(positionList) == 0:
                    raise SystemError("No overlaps found between {}s restricted by 'overlap'".format(btype))
            ## get new exon list with real genomic coordinates
            newBlockList = []
            prev = min(positionList) - 1
            blockStart = prev + 1
            for pos in positionList:
                if pos == positionList[-1]:
                    ## if the loop reach the end, store the last block
                    if pos - prev > 1:
                        blockEnd = prev
                        newBlockList.append([blockStart, blockEnd])
                        newBlockList.append([pos, pos])
                    else:
                        blockEnd = pos
                        newBlockList.append([blockStart, blockEnd])
                else:
                    if pos - prev > 1:
                        ## if the position is not continuous with previous, store the current block and start a new block
                        blockEnd = prev
                        newBlockList.append([blockStart, blockEnd])
                        blockStart = pos
                        prev = pos
                    else:
                        prev = pos
            ## remove the single number in newBlockList, like the "8" in [1,2,3,4,8,10,11]
            ## newBlockList: [[1,4], [8,10]]
            newBlockList = list(filter(lambda x:x[1] - x[0] > 0, newBlockList))
            if self.__overlap and btype == 'exon':
                if len(newBlockList) == 0:
                    raise SystemError("No overlaps found between blocks restricted by 'overlap'")
        else:
            ## if tx mode is on, should be very careful, do not destroy the trascript structure
            ## make intersected block in b as block-b
            ## if block-b overlap with the internal block of a, that they must be the same
            ## if the block-i overlap with the fisrt or last block of a, then the block-i should not overhang in the
            ## internal of the block of a
            b1len = len(block1List)
            b2len = len(block2List)
            newBlockList = []
            ## if block1 and block2 both have only 1 block
            if b1len == 1 and b2len == 1:
                bed1Row = ['chr1', block1List[0][0], block1List[0][1], 'bed1', 0, '+']
                bed2Row = ['chr1', block2List[0][0], block2List[0][1], 'bed2', 0, '+']
                bed1 = bed6ops(bed1Row)
                if self.__fun == 'merge':
                    bedops = bed1.merge(bed1Row, d=-1)
                elif self.__fun == 'intersect':
                    bedops = bed1.intersect(bed1Row)
                if bool(bedops) is True:
                    newBlockList.append([bedops.a.start, bedops.a.end])
                else:
                    raise SystemError("No {} (no overlap) when tx mode is on".format(self.__fun))
            elif b1len == 1 or b2len == 1:
                ## one of inputs has only 1 block
                ## always let block2List be blocksize == 1
                if b1len == 1:
                    block2List, block1List = block1List, block2List
                    b1len = len(block1List)
                    b2len = len(block2List)
                bed2Row = ['chr1', block2List[0][0], block2List[0][1], 'bed2', 0, '+']
                bed2 = bed6ops(bed2Row)
                overlapCount = 0
                overlapIndex = -1
                for i in range(b1len):
                    bed1Row = ['chr1', block1List[i][0], block1List[i][1], 'bed1', 0, '+']
                    if self.__fun == 'merge':
                        bedopsTmp = bed2.merge(bed1Row, d=-1)
                    elif self.__fun == 'intersect':
                        bedopsTmp = bed2.intersect(bed1Row)
                    if bool(bedopsTmp) is True:
                        overlapIndex = i
                        bedops = bedopsTmp
                        overlapCount += 1
                if overlapIndex < 0:
                    raise SystemError("No {} (no overlap) when tx mode is on".format(self.__fun))
                elif overlapCount > 1:
                    raise SystemError("2 or more overlaps found when tx mode is on".format(self.__fun))
                else:
                    if (btype == 'cds' or self.__part is True) and self.__fun == 'intersect':
                        newBlockList = [[bedops.a.start, bedops.a.end]]
                    else:
                        newBlockList = block1List
                        if overlapIndex == 0:
                            ## overlap with first block of a
                            if self.__part is False:
                                if block1List[0][1] != block2List[0][1]:
                                    raise SystemError("block right edges should be the same when tx mode is on".format(self.__fun))
                            else:
                                if self.__fun == 'intersect':
                                    newBlockList = [[bedops.a.start, bedops.a.end]]
                                else:
                                    newBlockList[0] = [bedops.a.start, bedops.a.end]
                        elif overlapIndex == (b1len - 1):
                            ## overlap with last block of a
                            if self.__part is False:
                                if block1List[-1][0] != block2List[0][0]:
                                    raise SystemError("block left edges should be the same when tx mode is on".format(self.__fun))
                            else:
                                if self.__fun == 'intersect':
                                    newBlockList = [[bedops.a.start, bedops.a.end]]
                                else:
                                    newBlockList[-1] = [bedops.a.start, bedops.a.end]
                        else:
                            ## ignore the internal overlap block to keep the transcript structure
                            if block1List[overlapIndex][0] <= block2List[0][0] and block2List[0][1] <= block1List[overlapIndex][1]:
                                if self.__part is False:
                                    raise SystemError("internal overlap is not allow({}) with part=False when tx mode is on".format(self.__fun))
                            else:
                                raise SystemError("internal overlap is not allow({}) when tx mode is on".format(self.__fun))
            else:
                ## a and b both have at least 2 blocks
                ## always make block2 has the smallest number of blocks
                switchFlag = False
                if b2len > b1len:
                    block2List, block1List = block1List, block2List
                    b1len = len(block1List)
                    b2len = len(block2List)
                    switchFlag = True
                ## record intersected block pairs
                block1Dict = {}
                block2Dict = {}
                for i in range(b1len):
                    bed1Row = ['chr1', block1List[i][0], block1List[i][1], 'bed1', 0, '+']
                    bed1 = bed6ops(bed1Row)
                    for j in range(b2len):
                        bed2Row = ['chr1', block2List[j][0], block2List[j][1], 'bed2', 0, '+']
                        if self.__fun == 'merge':
                            bedops = bed1.merge(bed2Row, d=-1)
                        elif self.__fun == 'intersect':
                            bedops = bed1.intersect(bed2Row)
                        if bool(bedops) is True:
                            if i not in block1Dict:
                                block1Dict[i] = {j:bedops}
                            else:
                                raise SystemError("2 or more overlaps on same block found when tx mode is on".format(self.__fun))
                            if j not in block2Dict:
                                block2Dict[j] = {i:bedops}
                            else:
                                raise SystemError("2 or more overlaps on same block found when tx mode is on".format(self.__fun))
                ## decode intersected blocks
                oblock1IndexList = sorted(block1Dict.keys())
                oblock2IndexList = sorted(block2Dict.keys())
                if len(oblock1IndexList) == 0:
                    raise SystemError("No {} (no overlap) when tx mode is on".format(self.__fun))
                firstA = oblock1IndexList[0]
                firstB = oblock2IndexList[0]
                lastA = oblock1IndexList[-1]
                lastB = oblock2IndexList[-1]
                ## check if intersect blocks are continuous
                block1ConFlag = (oblock1IndexList[0] + len(oblock1IndexList) - 1 == oblock1IndexList[-1])
                block2ConFlag = (oblock2IndexList[0] + len(oblock2IndexList) - 1 == oblock2IndexList[-1])
                if block1ConFlag is False or block2ConFlag is False:
                    raise SystemError("The overlaps is not continuous when tx mode is on".format(self.__fun))
                ## when only internal overlaps
                if self.__part is False or self.__fun == 'merge':
                    ## if part is False when running intersect(), or when running merge()
                    ## the overlaps should be started at the first block of a or b
                    ## ther overlaps should be ended at the last block of a or b
                    if ( firstA != 0 and firstB != 0 ) or ( lastA != (b1len - 1) and lastB != (b2len - 1) ):
                        raise SystemError("only internal overlaps is not allow when tx mode is on".format(self.__fun))
                ## fisrt overlap block
                if block1List[firstA][1] != block2List[firstB][1]:
                    if self.__rescue is False:
                        raise SystemError("block right edges should be the same when tx mode is on".format(self.__fun))
                    else:
                        bedops = block1Dict[firstA][firstB]
                        if switchFlag is True:
                            newBlockList.append([bedops.a.start, block1List[firstA][1]])
                        else:
                            newBlockList.append([bedops.a.start, block2List[firstB][1]])
                else:
                    bedops = block1Dict[firstA][firstB]
                    newBlockList.append([bedops.a.start, bedops.a.end])
                ## last overlap block
                if firstA < lastA:
                    if block1List[lastA][0] != block2List[lastB][0]:
                        if self.__rescue is False:
                            raise SystemError("block left edges should be the same when tx mode is on".format(self.__fun))
                        else:
                            bedops = block1Dict[lastA][lastB]
                            if switchFlag is True:
                                newBlockList.append([block1List[lastA][0], bedops.a.end])
                            else:
                                newBlockList.append([block2List[lastB][0], bedops.a.end])
                    else:
                        bedops = block1Dict[lastA][lastB]
                        newBlockList.append([bedops.a.start, bedops.a.end])
                ## check internal overlap blocks
                if lastA > firstA + 1:
                    for i in range(firstA + 1, lastA):
                        j = sorted(block1Dict[i].keys())[0]
                        if block1List[i][0] != block2List[j][0] or block1List[i][1] != block2List[j][1]:
                            if self.__rescue is False:
                                raise SystemError("internal blocks should be the same when tx mode is on".format(self.__fun))
                            else:
                                bedops = block1Dict[i][j]
                                if switchFlag is True:
                                    newBlockList.append(block1List[i])
                                else:
                                    newBlockList.append(block2List[j])
                        else:
                            newBlockList.append(block1List[i])
                ## for merge peaks
                if self.__fun == 'merge':
                    remainA = list(range(0, firstA)) + list(range(lastA + 1, b1len))
                    for i in remainA:
                        newBlockList.append(block1List[i])
                    remainB = list(range(0, firstB)) + list(range(lastB + 1, b2len))
                    for i in remainB:
                        newBlockList.append(block2List[i])
        ## get blockCount, blockSizes, blockStarts
        newBlockList = sorted(newBlockList, key=lambda x:x[0])
        blockCount = len(newBlockList)
        newStart = min(map(lambda x:x[0], newBlockList))
        newEnd = max(map(lambda x:x[1], newBlockList))
        blockSizesList = []
        blockStartsList = []
        for exon in newBlockList:
            blockSize = exon[1] - exon[0]
            blockStart = exon[0] - newStart
            blockSizesList.append(blockSize)
            blockStartsList.append(blockStart)
        ## joint the block size and block starts with ','
        blockSize = ','.join(map(str, blockSizesList))
        blockStarts = ','.join(map(str, blockStartsList))
        ## return final results
        return [newStart, newEnd, blockCount, blockSize, blockStarts]
    # calculated score
    def __getScore(self, scoreA, scoreB):
        socreList = [scoreA, scoreB]
        if self.__score == 'sum':
            score = sum(socreList)
        elif self.__score == 'min':
            score = min(socreList)
        elif self.__score == 'max':
            score = max(socreList)
        elif self.__score == 'average':
            score = sum(socreList) / 2
        elif self.__score == 'a':
            score = self.a.score
        elif self.__score == 'b':
            score = self.b.score
        else:
            raise SystemError("socre should be one of [sum, min, max, average, a, b]")
        return score
    # merge 2 bed12 based on exons
    def merge(self, b, score='sum', s=False, tx=True, overlap=False, cds=True):
        ## if overlap is True, only merge block when any overlaps found
        ## return a bed12ops object
        self.strand = s
        self.__tx = tx
        self.__score = score
        self.__fun = 'merge'
        self.__overlap = overlap
        self.__part = True
        self.__rescue = False
        ## get structures for self.a and self.b: exon, intron, cds, utr5, utr3
        self.b = buildbed(b)
        self.b = self.b.decode()
        self.clear = self.__check()
        chrom = self.a.chr
        name = '|'.join([self.a.name, self.b.name])
        if self.a.strand == self.b.strand:
            strand = self.a.strand
        else:
            strand = '.'
        ## how to get score
        newScore = self.__getScore(self.a.score, self.b.score)
        ## for exons
        try:
            start, end, bcount, bsize, bstart = self.__squeezeBlock(self.a.exon, self.b.exon, 'exon')
        except (SystemError, ValueError) as e:
            return False
        ## for cds
        if cds is True:
            try:
                tstart, tend, __, __, __ = self.__squeezeBlock(self.a.cds, self.b.cds, 'cds')
            except SystemError as e:
                tstart = start
                tend = start
        else:
            tstart = start
            tend = start
        ## get the final merged results
        row = [chrom, start, end, name, newScore, strand, tstart, tend, 255, bcount, bsize, bstart]
        self = bed12ops(row)
        return self
    # intersect 2 bed12 based on exons
    def intersect(self, b, score='sum', s=False, tx=True, part=False, cds=True, rescue=False):
        ## return a bed12ops object
        self.strand = s
        self.__tx = tx
        self.__score = score
        self.__fun = 'intersect'
        self.__overlap = True
        self.__part = part
        self.__rescue = rescue
        ## if part is True: for a has >2 blocks, b has >= 1 block, then return the overlaps instead of block-a or raise errors
        ## part only works with tx=True
        ## if rescue is True, then take b as a template to correct the blocks in a
        ## get structures for self.a and self.b: exon, intron, cds, utr5, utr3
        self.b = buildbed(b)
        self.b = self.b.decode()
        self.clear = self.__check()
        chrom = self.a.chr
        name = '|'.join([self.a.name, self.b.name])
        if self.a.strand == self.b.strand:
            strand = self.a.strand
        else:
            strand = '.'
        ## how to get score
        newScore = self.__getScore(self.a.score, self.b.score)
        try:
            start, end, bcount, bsize, bstart = self.__squeezeBlock(self.a.exon, self.b.exon, 'exon')
        except (SystemError, ValueError) as e:
            return False
        ## for cds
        if cds is True:
            try:
                tstart, tend, __, __, __ = self.__squeezeBlock(self.a.cds, self.b.cds, 'cds')
            except SystemError as e:
                tstart = start
                tend = start
        else:
            tstart = start
            tend = start
        ## get the final overlapd results
        row = [chrom, start, end, name, newScore, strand, tstart, tend, 255, bcount, bsize, bstart]
        try:
            self = bed12ops(row)
        except SystemError as e:
            return False
        return self
