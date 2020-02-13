"""
Constraints allow to find the density of the lowest state that satiffies certain conditions
on the distribution of charge over molecular fragments. If the charge distribution is chosen cleverly the constraint solution might even approximate an excited charge transfer (CT) state.
By enumerating all possible combinations of constraints (e.g. a charge sitting on any of the monomers in a conducting polymer such as polypyrrole ) we can get a coarse-grained set of potential energy surfaces (one state per constraint). Non-adiabatic dynamics within these states should be able
to simulate charge transport.
"""

class ChargeConstraint:
    """
    class holds the definition of the molecular fragments
    and the partial charges that are fixed to them.
    """
    def __init__(self, atomlist):
        """

        """
        self.atomlist = atomlist
        self.regions = []
    def addRegion(self, atomic_indeces, excess_charge):
        """
        Define a region by the indeces of the atoms that belong to it.
        Indeces refer to the atomlist with which the ChargeConstrain was initialized.
        """
        self.regions.append( (atomic_indeces, excess_charge) )
    def checkValidity(self):
        """
        check if all constraints can be fulfilled in principle or if there
        are conflicts.
        """
    def evaluateConstraints(self, dq):
        pass
    def gradientOfConstraint(self, dq):
        pass
