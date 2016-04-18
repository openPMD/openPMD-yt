"""
openPMD data structures


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from .fields import openPMDFieldInfo

from yt.utilities.file_handler import \
    HDF5FileHandler

import h5py
import numpy as np
import os
from yt.utilities.logger import ytLogger as mylog

class openPMDBasePathException(Exception) :
    pass

class openPMDBasePath :
    def  _setBasePath(self, handle) :
        """
        Set the base path for the first iteration found in the file.
        TODO implement into distinct methods:
            - __init__(self, handle)
            - getIterations(self)
            - getBasePath(self, iteration)
        """
        # basePath is fixed in openPMD 1.X to `/data/%T/`
        dataPath = u"/data"

        # if the file messed up the base path we avoid throwing a cluttered
        # exception below while looking for iterations:
        if handle.attrs["basePath"].decode("utf-8") != u"/data/%T/" :
            raise openPMDBasePathException("openPMD: basePath is non-standard!")

        # does `/data/` exist?
        if not u"/data" in handle :
            raise openPMDBasePathException("openPMD: group for basePath does not exist!")

        # find iterations in basePath
        list_iterations = []
        for i in list(handle[dataPath].keys()) :
            list_iterations.append(i)
        mylog.warning("openPMD: found {} iterations in file".format(len(list_iterations)))

        # We found no iterations in basePath
        # TODO in the future (see above) this can be a mylog.warning instead of an error
        if len(list_iterations) == 0 :
            raise openPMDBasePathException("openPMD: no iterations found in basePath!")

        # just handle the first iteration found
        mylog.warning("openPMD: only choose to load first iteration in file")
        self.basePath = "{}/{}/".format(dataPath, list_iterations[0])




class openPMDGrid(AMRGridPatch):
    """
    This class defines the characteristics of the grids
    Actually there is only one grid for the whole simolation box
    """
    _id_offset = 0
    __slots__ = ["_level_id"]

    def __init__(self, id, index, level=-1):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        # There is only one grid. So there are no parent or child grids
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "openPMDGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class openPMDHierarchy(GridIndex, openPMDBasePath):
    """
    Defines which fields and particles are created and read from the hard disk
    Furthermore it defines the characteristics of the grids
    """
    grid = openPMDGrid

    def __init__(self, ds, dataset_type='openPMD'):
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.dataset = ds
        self.index_filename = ds.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._setBasePath(self.dataset._handle)
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        """
        Test which fields and particle fields are in the file and defines the
        names of fields in yt

        This needs to set a self.field_list that contains all the available,
        on-disk fields.
        NOTE: Each should be a tuple, where the first element is the on-disk
        fluid type or particle type.  Convention suggests that the on-disk
        fluid type is usually the dataset_type and the on-disk particle type
        (for a single population of particles) is "io".
        look for fluid fields
        """
        meshesPath = self.dataset._handle.attrs["meshesPath"]
        particlesPath = self.dataset._handle.attrs["particlesPath"]
        output_fields = []

        # Read the names of the non-particle fields and add them to output_fields
        for group in self.dataset._handle[self.basePath + meshesPath].keys():
            try:
                for direction in self.dataset._handle[self.basePath + meshesPath + group].keys():
                    output_fields.append(group + "_" + direction)
            except:
                output_fields.append(group)

        # add the names to field_list
        self.field_list = [("openPMD", str(c)) for c in output_fields]

        # look for particle fields
        particle_fields = []

        # The following code read the particle type and the name of the particle
        # field out of the file
        if self.basePath + particlesPath in self.dataset._handle:
           for particle_type in self.dataset._handle[self.basePath + particlesPath].keys():
              for group in self.dataset._handle[self.basePath + particlesPath + particle_type].keys():

                 try:

                    key = self.dataset._handle[
                            self.basePath +
                            particlesPath +
                            particle_type +
                            "/" +
                            group].keys()
                    if key == []:
                       particle_fields.append(particle_type + "_" + group)
                       pass
                    else:
                       for direction in key:
                          particle_fields.append(
                             particle_type +
                             "_" +
                             group +
                             "_" +
                             direction)
                          pass

                 except:

                    particle_fields.append(particle_type + "_" + group)
                    # The name of the particle field and the particle type are added to
                    # field_list
                    self.field_list.extend([(str(c).split("_")[0], str(c).replace(str(c).split("_")[0], "particle"))
                                            for c in particle_fields])

    def _count_grids(self):
        """
        This needs to set self.num_grids
        Actually there is only one big grid
        """
        self.num_grids = 1

    def _parse_index(self):
        """
        The dimensions and the size of the grid is defined
        Actually there is only one grid so it has the same size like the
        simulationsbox
        """
        # This needs to fill the following arrays, where N is self.num_grids:
        meshesPath = self.dataset._handle.attrs["meshesPath"]
        particlesPath = self.dataset._handle.attrs["particlesPath"]

        self.grid_left_edge[
            0] = self.dataset.domain_left_edge  # (N, 3) <= float64
        self.grid_right_edge[
            0] = self.dataset.domain_right_edge  # (N, 3) <= float64
        self.grid_dimensions[
            0] = self.dataset.domain_dimensions  # (N, 3) <= int
# TODO this disables particle reads for now
#      Should might be set in _read_particles for each species,
#      also each species might need its own grid (?)
#        self.grid_particle_count[
#            0] = self.dataset._handle[
#                self.basePath +
#                particlesPath +
#                "/electrons/position/x"].shape[
#            0]  # (N, 1) <= int
        # self.grid_levels = 1           #(N, 1) <= int
        # self.grids = np.empty(1, dtype='object') #(N, 1) <= grid objects

        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype='object')
        # You have to inalize the grids
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels[i, 0])

    def _populate_grid_objects(self):
        """
        This function initializes the grids

        For each grid, this must call:
        grid._prepare_grid()
        grid._setup_dx()
        This must also set:
          grid.Children <= list of child grids
          grid.Parent   <= parent grid
        This is handled by the frontend because often the children must be
        identified.
        """

        # self._reconstruct_parent_child()

        for i in range(self.num_grids):
            self.grids[i]._prepare_grid()
            self.grids[i]._setup_dx()
        self.max_level = 0


class openPMDDataset(Dataset, openPMDBasePath):
    """
    A dataset object contains all the information of the simulation and
    is intialized with yt.load()
    
    TODO Ideally, a data set object should only contain a single data set.
         afaik, yt.load() can load multiple data sets and also supports
         multiple iteration-loading if done that way, e.g., from a prefix
         of files.
    """
    _index_class = openPMDHierarchy
    _field_info_class = openPMDFieldInfo

    def __init__(self, filename, dataset_type='openPMD',
                 storage_filename=None,
                 units_override=None):
        # This defines which of the data sets are meshes (fields)
        self.fluid_types += ('openPMD',)
        # This defines which of the data sets are particles
        self.particle_types = ["electrons", "ions", "all"]
        self.particle_types = tuple(self.particle_types)
        self.particle_types_raw = self.particle_types

        # Opens a HDF5 file and stores its file handle in _handle
        # All _handle objects refers to the file
        self._handle = HDF5FileHandler(filename)
        self._setBasePath(self._handle)
        Dataset.__init__(self, filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename


    def _set_code_unit_attributes(self):
        """
        This function defines the unit system of the on disk file
        """
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        self.length_unit = self.quan(1.0, "m")
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")

        #
        # These can also be set:
        self.velocity_unit = self.quan(1.0, "m/s")
        self.magnetic_unit = self.quan(1.0, "T")

    def _parse_parameter_file(self):
        """
        The parameters of simulation are loaded
        """
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be updated to be in code units at a later time.  This includes
        # the cosmological parameters.

        # read parameters out .h5 file
        f = self._handle

        meshesPath = f.attrs["meshesPath"]
        particlesPath = f.attrs["particlesPath"]
        positionPath = self.basePath + particlesPath + "/electrons/position/"

        # This defines the size of the simulaion box
        # TODO !!! The size is actually hardcoded
        self.unique_identifier = 0  # no identifier
        self.parameters = 0  # no additional parameters  <= full of code-specific items of use
        self.domain_left_edge = np.array(
            [-1.49645302e-05,
             -1.00407931e-06,
             -3.48883032e-06])  # <= array of float
        self.domain_right_edge = np.array(
            [1.49645302e-05,
             1.41565351e-06,
             1.54200006e-05])  # <= array of float64
        self.dimensionality = 3  # <= int

        fshape = []
        for i in range(3):
            try:
                fshape.append(f[self.basePath + meshesPath + "/B/x"].shape[i])
            except:
                fshape.append(1)
        self.domain_dimensions = np.array(
            [fshape[0],
             fshape[2],
                fshape[1]],
            dtype="int64")  # <= array of int64

        self.periodicity = (
            False,
            False,
            False)  # <= three-element tuple of booleans
        self.current_time = f[self.basePath].attrs[
            "time"]  # <= simulation time in code units
        self.refine_by = 2

        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.

        # Not a cosmological simulation
        self.cosmological_simulation = 0  # <= int, 0 or 1
        self.current_redshift = 0  # <= float
        self.omega_lambda = 0  # <= float
        self.omega_matter = 0  # <= float
        self.hubble_constant = 0  # <= float

    @classmethod
    def _is_valid(self, *args, **kwargs):
        """
        This function test if the (with yt.load()) a file could be opened with
        this frontend
        """
        # check if openPMD standard is used
        # return True
        # Test if the HDF5-file correspond to the openPMD standard
        need_attributes = ['openPMD', 'openPMDextension']

        valid = True
        try:
            fileh = h5py.File(args[0], mode='r')
            for na in need_attributes:
                if na not in fileh["/"].attrs.keys():
                    valid = False
            fileh.close()
        except:
            valid = False
        # if True is returned the file could be opened with the openPMD
        # standard
        return valid
