class RecordUniProt:

    def __init__(self, name, id, locale, status, fmol, bio, function2, leng, ec, seq = "-", evalu = "-", score = "-"):
        self.name = name
        self.id = id
        self.locale = locale
        self.status = status
        self.fmol = fmol
        self.bio = bio
        self.function2 = function2
        self.leng = leng
        self.ec = ec
        self.seq = seq
        self.evalu = evalu
        self.score = score



    def getAll(self):
        return self.name, self.id, self.locale, self.status, self.fmol, self.function2, self.leng, self.ec, self.seq, self.evalu, self.score

    def getName(self):
        return self.name

    def getId(self):
        return self.id

    def getlocale(self):
        return self.locale

    def getstatus(self):
        return self.status

    def getfmol(self):
        return self.fmol

    def getbio(self):
        return self.bio

    def getfunction2(self):
        return self.function2

    def getleng(self):
        return self.leng

    def getEC(self):
        return self.ec

    def getSeq(self):
        return self.seq

    def getEvalu(self):
        return self.evalu

    def getScore(self):
        return self.score
