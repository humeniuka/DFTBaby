"""
Read the dynamics.out files that are produced with Jens' program using the "Full" output option
"""

############ INPUT ##################

from DFTB import AtomicData            

class _Block(object):
    def __init__(self, args=[]):
        """called at the beginning of a block"""
        pass
    def process(self, line):
        pass
    def close(self):
        """called at the end of a block"""
        pass
    def getData(self):
        """returns a nice representation of the data read"""
        pass

class _IgnoreBlock(_Block):
    pass

class time_Block(_Block):
    def __init__(self, args=[]):
        self.atomlist = []
        self.velocities = []
        self.read = None
        self.time = float(args[0].strip()[:-2])
    def process(self, l):
        # Which kind of information should be read?
        l = l.strip()
        if self.read == None:
            if l == "coordinates":
                self.read = "C"
                return
        elif self.read == "C":
            if l == "velocities":
                self.read = "V"
                return
        # Read coordinates or velocities line by line
        if self.read == "C":
            atom,X,Y,Z = l.split()
            Zi = AtomicData.atomic_number(atom)
            # positions in bohr
            pos = float(X), float(Y), float(Z)
            self.atomlist.append( (Zi,pos) )
        elif self.read == "V":
            vx,vy,vz = l.split()
            vel = float(vx), float(vy), float(vz)
            self.velocities.append( vel )
    def close(self):
        pass
    def getData(self):
        return (self.time, self.atomlist, self.velocities)

class state_of_next_step(_Block):
    def __init__(self, args=[]):
        self.step = int(args[0])
    def getData(self):
        return self.step

def parse_dynamics_out(out_file):
    fh = open(out_file, "r")
    while True:
        l = fh.readline()
        if l == "":
            # end of file
            break
        block = _IgnoreBlock()
        for l in fh:
            if ":" in l: # 
                args = l.split(":")
                block_type = args[0].strip()
                # remove strange characters
                block_type = block_type.replace(" ", "_")
                # close previous block
                block.close()
                if block.__class__ != _IgnoreBlock:
                    yield block
                try:
                    block = eval(block_type + "_Block")(args[1:])
                except NameError as e:
#                    print e
#                    print "Block %s not implemented" % block_type 
                    block = _IgnoreBlock()
            else:
                block.process(l)
    block.close()
    if block.__class__ != _IgnoreBlock:
        yield block

def t_geometries_velocities(dynout):
    for b in parse_dynamics_out(dynout):
        if b.__class__ == time_Block:
            t,atomlist,velocities = b.getData()
            yield (t,atomlist,velocities)

################### OUTPUT ###################

def write_initial_conditions(atomlist, velocities, in_file):
    """
    write initial conditions
    """
    Nat = len(atomlist)
    txt = "%d\n" % Nat
    # positions in bohr
    for Zi, (x,y,z) in atomlist:
        atname = AtomicData.atom_names[Zi-1]
        txt += "%4s    %+10.15f    %+10.15f    %+10.15f\n" % (atname, x, y, z)
    # velocities in bohr/s
    for i in range(0, Nat):
        vx,vy,vz = velocities[i]
        txt += "   %+10.15f    %+10.15f   %+10.15f\n" % (vx,vy,vz)
    #
    fh = open(in_file, "w")
    print>>fh, txt

if __name__ == "__main__":
    import sys
    for b in parse_dynamics_out(sys.argv[1]):
        if b.__class__ == time_Block:
            t,atomlist,velocities = b.getData()
            print t


