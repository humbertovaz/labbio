class Record:

    def __init__(self,locustag,genename,strand):
        self.locustag = locustag
        self.genename = genename
        self.strand = strand

    def getlocustag(self):
        return self.locustag

    def getgenename(self):
        return self.genename

    def getstrand(self):
        return self.strand

    def getAll(self):
        return self.locustag, self.genename, self.strand




