Developer readme
updated 7/11/2017 by dve

Software structure:
All simulations require exactly one CList, MutationHandler, and SimParams. Each simulation also requires at least one CellType and at least one Clone.
-CList is main simulator container class, and represents a population of cells (or organisms). It encodes how cells will be chosen for reproduction and when they reproduce. It manages simulation time and the cell type space in the population. If you want to implement a new birth death process, you should start by subclassing a class from this family. The advance() method is the key method here; it is called repeatedly to run the simulation, and updates the CList to reflect changes that happened during that time step. CLists contain CellTypes, which in turn contain Clones.
-CellTypes are groups of cells that share a specific heritable characteristic (eg genotype, barcode, etc). When cells of a certain CellType reproduce, their progeny will be of the same CellType unless a mutation occurs. Clones are contained within each CellType, and different clones within a single CellType may have different birth rates. The CellType phylogenetic tree (which CellTypes are descendents of other CellTypes through mutation) is also recorded.
-Clones represent individual cells or groups of cells with the same CellType and same birth rate. The key method in this class is the reproduce() method, which defines what happens when a Clone of that class reproduces. The reproduce() method defines the inheritance model and birth rate distribution of the daughter cells. reproduce() also determines whether a mutation occurs during that division, but what happens if a mutation occurs is determined by the MutationHandler of the simulation.
-MutationHandlers specify the mutational behavior of the population. More generally, they define what transitions between CellTypes may occur and at what frequency. They are largely responsible for managing the type space and type-specific fitness differences. If a new CellType should be created after a mutation, the MutationHandler makes a new CellType, inserts it into the CList, and determines the birth rate of the new cell based on its new type.
-SimParams holds information about the simulation: how many trials to run, when to stop, and what data to collect. This class is also responsible for initializing instances of other classes (Clones, MutationHandlers, etc).

The following classes are not strictly required to run a simulation but add functionality.
-OutputWriters determine what data will be written to output files before, during, and after each trial of the simulation. The key methods are beginAction(), duringSimAction(), and finalAction(), which determine what (if anything) will be written to output files before, during, and after each trial.
-EndListeners are conditions (in addition to the no Clone condition) that end simulation trials. Popular examples include stopping a run after a certain number of time steps, or once a certain CellType appears in the population.

Developers will largely be adding extra CList, Clone, MutationHandler, OutputWriter, and EndListener classes. I strongly recommend leaving the rest of the architecture alone.

Known issues:
-should fix hierarchy (and name) of CList/MoranPop. Both a branching process simulator and a Moran simulator should inherit from a virtual population class.
-speed could be increased by storing clones in a b-tree like structure. you would be able to reduce reproduction time from O(n) to O(log n).
