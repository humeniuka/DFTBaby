# NOT FINISHED
"""
save results of a TD-DFTB calculation in a text format
reminiscent of the Gaussian checkpoint file. Each block
begins with a line like
  $forces float (10,3)
followed by the numerical data
"""

class _Block(object):
    def __init__(self):
        return self
    def write(self, txt, data):
        pass
    def read(self, txt):
        pass

class _numpy_Block(_Block):
    def write(self, txt, data, **kwds):
        sh = data.shape
        arr = data.ravel()
        dtype = str(arr.dtype)
        tag = self.__class__.__name__[1:].replace("_Block", "")
        txt += "$%s type=%s shape=%s" % (tag, dtype, sh.replace(" ",""))
        n = len(arr)
        for i in range(0, n):
            txt += "%s " % arr[i],
            if i % 5 == 0:
                txt += "\n"
        return txt
    def read(self, tag_line, txt):
        tag, dtype, shape = tagline[1:].strip().split()
        arr = map(float, txt.split())
        data = np.reshape(np.array(arr, dtype=dtype), tuple(shape))
        return tag, data

class _gradient_Block(_Block):
    pass

def get_next_block(fh):
    buf = ""
    cur_block = None
    while True:
        l = fh.readline()
        if l == "":
            # end of file
            if cur_block != None:
                yield cur_block.read(l, buf)
            raise StopIteration
        if l[0] == "$":
            # close previous block and start new block
            if cur_block != None:
                yield cur_block.read(l, buf)
                # reset buffer
                buf = ""
            # determine block type from tag name
            block_type = l.split()[0][1:]
            try:
                cur_block = eval(block_type + "_Block")()
            except:
                cur_block = None
        else:
            buf += l

def parseMyCheckpointFile(filename):
    fh = fh.open(filename)
    for (tag, data) in get_next_block(fh):
        yield tag, data
    fh.close()

