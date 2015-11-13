"""
OpenPMD-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.extern.six import u, b, iteritems
from contextlib import contextmanager
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.selection_routines import mask_fill, AlwaysSelector
import h5py
import numpy as np

_convert_mass = ("particle_mass","mass")


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# This class loads the data from the HDF5-file.
# TODO Data should be loaded chunk-wise to support parallelism. Fields can be chunked arbitrarily
#      in space, particles should be read via particlePatches.

class IOHandlerOpenPMD(BaseIOHandler):
    _dataset_type = "openPMD"
    _field_dtype = "float32"

    def __init__(self, ds, *args, **kwargs):

	self.ds = ds
    # ds._handle is a HDF5-file which is loaded with h5py.File(filename)
	self._handle = ds._handle
	self.basePath = self._handle["/"].attrs["basePath"]
	self.meshPath = self._handle["/"].attrs["meshesPath"]
	self.particlesPath = self._handle["/"].attrs["particlesPath"]

    # This function read the coords of the particles out of the file
    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
	#print "readParCoords"

	chunks = list(chunks)
	self._array_fields = {}
    # A chunk returns a grid element with a filename
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if g.filename is None: continue

                if f is None:
                    #print "Opening (count) %s" % g.filename

		    #open HDF5-file to read
                    f = h5py.File(g.filename, "r")

                dds = f.get(f.attrs["basePath"])
		for ptype, field_list in sorted(ptf.items()):
		    # ptf.items() returns a list with all known particle fields
            # and particle types

                    #print field_list
		    if field_list == "particle_mass":  # do not know if necessary
			continue
		    pds = dds.get("%s/%s" % (f.attrs["particlesPath"], ptype))


                    # read particle coords out of file
                    x, y, z = (np.asarray(pds.get("position/" + ax).value, dtype="=f8")
                               for ax in 'xyz')
                    for field in field_list:
            # Save the size of the particle fields
			nfield = field.replace("particle_","")
			nfield = nfield.replace("_","/")
			#print "Field " + nfield

                        if np.asarray(pds[nfield]).ndim > 1:
                            self._array_fields[field] = pds[nfield].shape
            # Returns the coords of the particle
                    yield ptype, (x, y, z)
            if f: f.close()

    # This function read the particle fields out of file
    def _read_particle_fields(self, chunks, ptf, selector):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
	#print "readPartField"
	chunks = list(chunks)
    # A chunk returns a grid element wiht a file name
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if g.filename is None: continue
                if f is None:
                    #print "Opening (count) %s" % g.filename

		    #open a HDF5-file to read
                    f = h5py.File(u(g.filename), 'r')
                    #print g.filename + " io.py"
                ds = f.get(f.attrs["basePath"])
		for ptype, field_list in sorted(ptf.items()):
            # ptf.items() returns a list of all known particle fields and
            # particle types
                    pds = ds.get("%s/%s/" % (f.attrs["particlesPath"],ptype))

            # The particle coords have to be loaded again
                    x, y, z = (np.asarray(pds.get("position/" + ax).value, dtype="=f8")
                               for ax in 'xyz')

                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None: continue
                    for field in field_list:
			nfield = field.replace("particle_","")
			nfield = nfield.replace("_","/")
            # Load the field informations out of file
            # !!! Incomplete because mass and charge are a attributes and
            # the rest are arrays
			if field == "particle_mass" or field == "particle_charge":
			    data = np.full(x.shape[0],pds.get(nfield).attrs["value"], "=f8")
			else:
			    data = np.asarray(pds.get(nfield), "=f8")
            # Here you could multiply mass with weighting
                        #if field in _convert_mass:
                        #    data *= g.dds.prod(dtype="f8")

            # This returns particle type, field name and the masked field data
                        yield (ptype, field), data[mask]
            if f: f.close()

    # This function reads the rest of the fields from the file
    def _read_fluid_selection(self, chunks, selector, fields, size):
        # This needs to allocate a set of arrays inside a dictionary, where the
        # keys are the (ftype, fname) tuples and the values are arrays that
        # have been masked using whatever selector method is appropriate.  The
        # dict gets returned at the end and it should be flat, with selected
        # data.  Note that if you're reading grid data, you might need to
        # special-case a grid selector object.
	rv = {}
        # Now we have to do something unpleasant
        #print "readFluid"
	chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
        # Open the file
            f = h5py.File(u(g.filename), 'r')
	    for ftype, fname in fields:
        # Data is loaded with the _read_data function
                rv[ftype, fname] = self._read_data(g, fname)

	    f.close()
	    #print rv
            return rv

        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))
        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
	#print "Reading %s cells of %s fields in %s grids",size, [f2 for f1, f2 in fields], ng

        ind = 0
        for chunk in chunks:
            for g in chunk.objs:
        # Open file
                f = h5py.File(u(g.filename), 'r')
	        for ftype, fname in fields:
        # Data is loaded with the _read_data function
		    rv[ftype, fname] = self._read_data(g, fname)

	        f.close()

        return rv

    # This function read the field data out of file
    def _read_data(self, grid, field):

	data = []
	if field.startswith("particle") :
        # particles are not loaded here
	    #data = self._handle[self.basePath+self.particlesPath+field.replace("_","/")]
	    pass
	else:
        # This reads the file
	    data = self._handle[self.basePath+self.meshPath+field.replace("_","/")]
        return np.array(data).flatten()


    # TODO This function is for caching. Actually this function does not work
    # For parallelism and working with big files it is recommended to use it
    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk, and is only used for
        # caching.
	#print "readChunk"
	f = self._handle
        rv = {}
        for g in chunk.objs:
            rv[g.id] = {}
        # Split into particles and non-particles
        fluid_fields, particle_fields = [], []
        for ftype, fname in fields:
            if ftype in self.ds.particle_types:
                particle_fields.append((ftype, fname))
            else:
                fluid_fields.append((ftype, fname))
        if len(particle_fields) > 0:
            selector = AlwaysSelector(self.ds)
            rv.update(self._read_particle_selection(
                [chunk], selector, particle_fields))
        if len(fluid_fields) == 0: return rv
        for field in fluid_fields:
            ftype, fname = field
            ds = f["/%s" % fname]
            ind = 0
            for gs in grid_sequences(chunk.objs):
                start = gs[0].id - gs[0]._id_offset
                end = gs[-1].id - gs[-1]._id_offset + 1
                data = ds[start:end,:,:,:].transpose()
                for i, g in enumerate(gs):
                    rv[g.id][field] = np.asarray(data[...,i], "=f8")
        return rv
