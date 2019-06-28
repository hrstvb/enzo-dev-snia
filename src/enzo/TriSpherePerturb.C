#include <stdio.h>
#include <math.h>
#include "DebugMacros.h"
#include "myenzoutils.h"

struct TriSphere
{
	size_t nV, nE, nF; // Number of vertices, edges, faces. V - E + F = 2
	double R, A, center[3], topBaseSize, bottomBaseSize; // radius, perturbation amplitude, center
	double *vx, *vy, *vz; // vertices coordinates
	size_t *ev1, *ev2; // edge vertices indices
	size_t *fe1, *fe2, *fe3; // face edges indices
	size_t nRefinements, nTimesRefined; // number of refinements
	bool hasBottomBase, hasTopBase, useTopAsOrigin, isHourglass;

	size_t cache_faceIndex;
	static const size_t cache_vectors_count = 24;

	double cache_vectors[cache_vectors_count][3];
	double (*_rect)[3];
	double (*_faceOrigins)[3];
	double (*_cofactors)[3];
	double (* const test) = cache_vectors[0]; //   3    0

	union
	{
		struct
		{
			double faceVs[3][3];
			double basesVs[6][3];
			double botCenter[3];
			double topCenter[3];
			double insidePoint[3];
			double rectL[3];
			double rectR[3];
			double origins[5][3];
			double cofactors[5][3];
		};
		double vector_cache[cache_vectors_count][3];
	};

	static void refineCounts(size_t* nV, size_t* nE, size_t* nF)
	{
		if(*nV <= 0)
		{
			// Icosahedron counts
			*nV = 12; //Vertices
			*nE = 30; //Edges
			*nF = 20; //Faces
		}
		else
		{
			// A one time refinement means that we divide each edge in two
			// and each face in four.
			*nV += *nE; // *e* new vertices for the middle of each edge.
			// Each edge is divided in two, which doubles the number of edges.
			// In addition three new edges appear on each face, connecting
			// the middles of the old edges.
			*nE = 2 * *nE + 3 * *nF;
			// divide each face in four, which quadruples the number of faces.
			*nF *= 4;
		}
	}

	/*
	 * Returns a TriSphere, starting from an icosahedron, refined as many times as
	 * the specified resolution, *dx* and *quality* permit.  The resulting triangular
	 * faces should have the smallest distance from the triangle center to any of the
	 * sides rectified by at least *quality* grid zones.  For an equilateral triangle
	 * with sides l, this means sqrt(3)/6 * l >= dx * quality.
	 */
	TriSphere(double radius, double amplitude, double dx, double quality, double bottomBaseSize, double topBaseSize)
	{
		double l0 = sqrt(1 + M_GOLDEN_RATIO * M_GOLDEN_RATIO) / 2; // Initial edge length for radius 1.
		l0 *= sqrt(3) / 6;
		l0 *= radius;
		double l = l0;
		l /= quality;

		nRefinements = 0;
		while((l /= 2) > dx)
			nRefinements++;

		TRACEF("l0=%e R=%e dx=%e q=%e l=%e -> %lld refinements", l0, radius, dx, quality, l, nRefinements);
		init(NULL, radius, amplitude, nRefinements, bottomBaseSize, topBaseSize);
	}

	TriSphere(double radius, double amplitude, size_t nRefinements, double bottomBaseSize, double topBaseSize)
	{
		init(NULL, radius, amplitude, nRefinements, bottomBaseSize, topBaseSize);
	}

	void init(double* center, double radius, double amplitude, int nRefinements, double bottomBaseSize,
		double topBaseSize)
	{
		vx = vy = vz = NULL;
		ev1 = ev2 = fe1 = fe2 = fe3 = NULL;
		nV = nE = nF = nTimesRefined = 0;
		R = radius;
		A = amplitude;
		this->topBaseSize = topBaseSize;
		this->bottomBaseSize = bottomBaseSize;
		if(center)
			arr_cpy(this->center, center, 3);
		else
			arr_set(this->center, 3, 0);

		hasBottomBase = bottomBaseSize != 0;
		hasTopBase = topBaseSize != 0;
		isHourglass = bottomBaseSize * topBaseSize < 0;
		useTopAsOrigin = fabs(topBaseSize) > fabs(bottomBaseSize);

		cache_faceIndex = -1;
		_rect = cache_vectors + 0;
		_cofactors = cache_vectors + 2;
		_faceOrigins = cache_vectors + 7;

//		arr_set(cache_vectors, 3 * cache_vectors_count, 0);

		cacheClear();

		this->nRefinements = nRefinements;
		for(int i = 0; i <= nRefinements; i++)
			refineCounts(&nV, &nE, &nF);

		allocateArrays(nV, nE, nF);
		initIcosahedron();
		for(int i = 1; i <= nRefinements; i++)
			refine();

		TRACEF("TriSphere: V,E,F= %lld %lld %lld after %lld refinements", nV, nE, nF, nRefinements);
	}

	~TriSphere()
	{
		freeArrays();
	}

	void allocateArrays(size_t nv, size_t ne, size_t nf)
	{
		vx = new double[nv];
		vy = new double[nv];
		vz = new double[nv];
		ev1 = new size_t[ne];
		ev2 = new size_t[ne];
		fe1 = new size_t[nf];
		fe2 = new size_t[nf];
		fe3 = new size_t[nf];
	}

	void freeArrays()
	{
		delete vx;
		delete vy;
		delete vz;
		delete ev1;
		delete ev2;
		delete fe1;
		delete fe2;
		delete fe3;

		vx = vy = vz = NULL;
		ev1 = ev2 = fe1 = fe2 = fe3 = NULL;
	}

	double getVertexRadius(size_t i)
	{
		return sqrt(square(vx[i]) + square(vy[i]) + square(vz[i]));
	}

	void setVertexRadius(size_t i, double r)
	{
		r /= getVertexRadius(i);
		vx[i] *= r;
		vy[i] *= r;
		vz[i] *= r;
	}

	/*
	 * Adds a vertex and returns its index.
	 * The radius-vector (*x*,*y*,*z*) is rescaled to *radius*.
	 */
	size_t addVertex(double x, double y, double z)
	{
		vx[nV] = x;
		vy[nV] = y;
		vz[nV] = z;
		setVertexRadius(nV, R);
		return nV++;
	}

	/*
	 * Gets the vertices' indices of an edge.
	 */
	void getEdge(size_t *i1, size_t *i2, size_t j)
	{
		*i1 = ev1[j];
		*i2 = ev2[j];
	}

	/*
	 * Sets the j-th edge and returns its index, j.
	 * May switch the two vertices so that i1 < i2.
	 */
	size_t setEdge(size_t j, size_t i1, size_t i2)
	{
		if(i1 < i2)
		{
			ev1[j] = i1;
			ev2[j] = i2;
		}
		else
		{
			ev1[j] = i2;
			ev2[j] = i1;
		}
		return j;
	}

	/*
	 * Adds and new edge and returns its index.
	 * May switch the two vertices so that i1 < i2.
	 */
	size_t addEdge(size_t i1, size_t i2)
	{
		return setEdge(nE++, i1, i2);
	}

	/*
	 * Gets the k-th face.
	 */
	void getFace(size_t* j1, size_t* j2, size_t* j3, size_t k)
	{
		*j1 = fe1[k];
		*j2 = fe2[k];
		*j3 = fe3[k];
	}

	/*
	 * Gets the k-th face.
	 */
	void getFace(double* v1, double* v2, double* v3, size_t k)
	{
		// if the k-th face's edges are j1(i1, i2), j2(i2, i3), j3(i1, i3), then
		// i1 < i2 < i3
		size_t j = fe1[k];
		size_t i = ev1[j];
		v1[0] = vx[i];
		v1[1] = vy[i];
		v1[2] = vz[i];
		j = fe2[k];
		i = ev1[j];
		v2[0] = vx[i];
		v2[1] = vy[i];
		v2[2] = vz[i];
		i = ev2[j];
		v3[0] = vx[i];
		v3[1] = vy[i];
		v3[2] = vz[i];
	}

	/*
	 * Sets the edges of the k-th face and returns its index, k>=0.
	 * May store the edge indices in a different order, so that
	 * if i1 < i2 < i3 are the three vertices, then the edges are
	 * fe1 = (i1,I2), fe2 = (i2,i3), fe3 = (i1,i3)
	 * Returns -1 if j1 and j2 are not connected.
	 * Returns -2 if the three edges don't make a triangle.
	 //# ...................
	 //#.........3..........
	 //#........./\.........
	 //#......../..\........
	 //#......fe3..fe2......
	 //#....../......\......
	 //#...../___fe1__\.....
	 //#....1..........2....
	 //#....................
	 */
	size_t setFace(size_t k, size_t j1, size_t j2, size_t j3)
	{
		size_t a, b, i31, i32;

		// Find if j1 and j2 have a common vertex.
		// Get the free ends in a and b.
		if(sortAngleEnds(&a, &b, j1, j2))
			return -1;

		i31 = ev1[j3];
		i32 = ev2[j3];
		// Does j3 close the triangle?
		if(i31 != a || i32 != b)
			return -2;

		sortEdges(&j1, &j2, &j3);
		fe1[k] = j1;
		fe2[k] = j3;
		fe3[k] = j2;
		return k;
	}

	/*
	 * Given three edge indices, adds a new face and returns its index.
	 * May store the edge indices in a different order. See setFace.
	 */
	size_t addFace(size_t j1, size_t j2, size_t j3)
	{
		int retval = setFace(nF, j1, j2, j3);
//		if(nTimesRefined>0)
//		TRACEF0("%lld %lld %lld %lld", retval, j1, j2, j3);
		if(retval < 0)
			return retval;
		return nF++;
	}

	double verticesDistance(size_t i1, size_t i2)
	{
		return sqrt(square(vx[i2] - vx[i1]) + square(vy[i2] - vy[i1]) + square(vz[i2] - vz[i1]));
	}

	double getEdgeLength(size_t j)
	{
		return verticesDistance(ev1[j], ev2[j]);
	}

	/*
	 * Returns -1 if edge j1 comes "before" j2, or returns +1
	 * if j1 comes "after" j2. Edges are comapred by their first
	 * vertex indices, and then - by the second.
	 */
	int cmpEdges(size_t j1, size_t j2)
	{
		int s = cmp(ev1[j2], ev1[j1]);
		if(s)
			return s;
		return cmp(ev2[j2], ev2[j1]);
	}

	/*
	 * Initializes the data with the vertices, edges and faces
	 * of an icosahedron inscribed in a sphere of given *radius*.
	 */
	void initIcosahedron()
	{
		const double a = M_GOLDEN_RATIO;
		nV = nE = nF = 0;

		// Initialize the vertices.
		// The first two must belong to the same edge.
		addVertex(0, 1, a);
		addVertex(0, -1, a);
		addVertex(0, 1, -a);
		addVertex(0, -1, -a);
		addVertex(a, 0, 1);
		addVertex(a, 0, -1);
		addVertex(-a, 0, 1);
		addVertex(-a, 0, -1);
		addVertex(1, a, 0);
		addVertex(-1, a, 0);
		addVertex(1, -a, 0);
		addVertex(-1, -a, 0);

		for(size_t i = 0; i < nV; i++)
		{
			vx[i] += center[0];
			vy[i] += center[1];
			vz[i] += center[2];
		}

		initEdgesByLength();
		initFaces();
	}

	/*
	 * Identify edges by length.
	 * Expects vertices 0 and 1 to belong to the same edge.
	 */
	void initEdgesByLength()
	{
		const double L = verticesDistance(0, 1);
		const double TOLLERANCE = 1e-6;

		for(size_t i1 = 0; i1 < nV - 1; i1++)
			for(size_t i2 = i1 + 1; i2 < nV; i2++)
			{
				if(fabs(L - verticesDistance(i1, i2)) < TOLLERANCE)
					addEdge(i1, i2);
			}
	}

	/*
	 * Returns 0 if the edges j1 and j2 have a common vertex and
	 * sets i1 and i2 to the other two "free" vertex indices,
	 * so that i1 < i2.
	 * Returns -1 if j1 and j2 are not connected.
	 */
	int sortAngleEnds(size_t* const i1, size_t* const i2, const size_t j1, const size_t j2)
	{
		// Get the four vertices.
		size_t i11 = ev1[j1];
		size_t i21 = ev1[j2];
		size_t i12 = ev2[j1];
		size_t i22 = ev2[j2];

		// Find the two common vertices.
		// The other two will end up in i12 and i22.
		if(i11 == i22)
		{
			i22 = i21;
		}
		else if(i12 == i21)
		{
			i12 = i11;
		}
		else if(i12 == i22)
		{
			i12 = i11;
			i22 = i21;
		}
		else if(i11 != i21)
			return -1; // N common vertice (not an angle)

		// sort the two indices.
		if(i12 < i22)
		{
			*i1 = i12;
			*i2 = i22;
		}
		else
		{
			*i1 = i22;
			*i2 = i12;
		}

		return 0;
	}

	/*
	 * Sorts edge indices so that, j1 <= j2 <= j3.
	 * Edges are comapred by their first vertex indices,
	 * then by the second.
	 */
	void sortEdges(size_t* const j1, size_t* const j2, size_t* const j3)
	{
		if(cmpEdges(*j1, *j2) < 0)
		{
			size_t temp = *j1;
			*j1 = *j2;
			*j2 = temp;
		}
		if(cmpEdges(*j2, *j3) < 0)
		{
			size_t temp = *j2;
			*j2 = *j3;
			*j3 = temp;
		}
		if(cmpEdges(*j1, *j2) < 0)
		{
			size_t temp = *j1;
			*j1 = *j2;
			*j2 = temp;
		}
	}

	void initFaces()
	{
		// Find the faces
		const size_t nEOld = nE;
		for(size_t j1 = 0; j1 < nEOld - 2; j1++)
			for(size_t j2 = j1 + 1; j2 < nEOld - 1; j2++)
				for(size_t j3 = j2 + 1; j3 < nEOld; j3++)
					addFace(j1, j2, j3);
	}

	/*
	 * Bisects the j-th edge:
	 * (i) creates a new vertex in the middle of the arc.
	 * (ii) replaces the old edge with first "half" edge.
	 * (iii) The second "half" is created at the E+j position.
	 */
	void bisectEdge(size_t j)
	{
		size_t i1 = ev1[j];
		size_t i2 = ev2[j];
		size_t i3 = addVertex((vx[i1] + vx[i2]) / 2, (vy[i1] + vy[i2]) / 2, (vz[i1] + vz[i2]) / 2);
		setEdge(j, i1, i3);
		addEdge(i2, i3);
	}

	/*
	 * Divides the k-th face in four new faces.
	 * Edges need to be bisected already.
	 * One of those goes in the k-th position, replacing the old face.
	 * The other three are appended.
	 * Three new edges are also created.
	 */
	void divideFace(const size_t k, const size_t nEOld)
	{
		//# Vertices, edges and faces during face refinement:
		//#..............................................
		//#......................3.......................
		//#....................../\......................
		//#...................../..\.....................
		//#..................../....\....................
		//#.................../......\...................
		//#.................36........35.................
		//#................./....k3....\.................
		//#................/............\................
		//#.............../..............\...............
		//#............6./_______56_______\.5............
		//#............./\................/\.............
		//#............/..\............../..\............
		//#.........../....\............/....\...........
		//#........../......\....k4..../......\..........
		//#........16.......46........45.......25........
		//#......../....k1....\....../....k2....\........
		//#......./............\..../............\.......
		//#....../..............\../..............\......
		//#...../_______14_______\/_______24_______\.....
		//#...1...................4..................2...
		//#..............................................
		const size_t j14 = fe1[k];
		const size_t j25 = fe2[k];
		const size_t j16 = fe3[k];
		const size_t j24 = nEOld + j14;
		const size_t j35 = nEOld + j25;
		const size_t j36 = nEOld + j16;
		//const size_t i1 = ev1[j14];
		const size_t i4 = ev2[j14];
		//const size_t i2 = ev1[j25];
		const size_t i5 = ev2[j25];
		//const size_t i3 = ev1[j36];
		const size_t i6 = ev2[j36];
		const size_t j45 = addEdge(i4, i5);
		const size_t j56 = addEdge(i5, i6);
		const size_t j46 = addEdge(i4, i6);

		/* size_t k1 = */setFace(k, j14, j46, j16);
		/* size_t k2 = */addFace(j24, j45, j25);
		/* size_t k3 = */addFace(j35, j56, j36);
		/* size_t k4 = */addFace(j45, j56, j46);
	}

	void refine()
	{
		size_t nEOld = nE;
		size_t nFOld = nF;

		for(size_t j = 0; j < nEOld; j++)
			bisectEdge(j);

		for(size_t k = 0; k < nFOld; k++)
			divideFace(k, nEOld);

		nTimesRefined++;
	}

	/*
	 * Returns the four vertices and the center of the tetrahedron.
	 */
	void getTetrahedron(size_t k, double v[10][3])
	{
		double *cb = (double*) (v + 7); // bottom center
		double *ct = (double*) (v + 8); // top center (spike tip)
		double *cc = (double*) (v + 6); // body (spike) center
		arr_set((double*) (v[0]), 30, 0);

		// if the k-th face's edges are j1(i1, i2), j2(i2, i3), j3(i1, i3), then
		// i1 < i2 < i3
		size_t j = this->fe1[k];
		size_t i = ev1[j];
		v[0][0] = vx[i];
		v[0][1] = vy[i];
		v[0][2] = vz[i];
		j = this->fe2[k];
		i = ev1[j];
		v[1][0] = vx[i];
		v[1][1] = vy[i];
		v[1][2] = vz[i];
		i = ev2[j];
		v[2][0] = vx[i];
		v[2][1] = vy[i];
		v[2][2] = vz[i];

		arr_set(ct, 6, 0);
		for(i = 0; i < 3; i++)
			arr_xpy(ct, v[i], 3);

		arr_axpy(cb, ct, 3, 1. / 3.); // base center
		scaleto(ct, 3, R + A); // spike tip
		arr_axpy(cc, cb, 3, 0.5);
		arr_axpy(cc, ct, 3, 0.5); // an inside point

		for(i = 0; i < 3; i++)
		{
			double a[3];
			arr_cpy(a, v[i], 3);
			arr_axpy(a, cb, 3, -1);
			arr_cpy(v[i + 3], ct, 3);
			if(topBaseSize > 0)
				arr_axpy(v[i + 3], a, 3, topBaseSize);
			if(bottomBaseSize > 0)
				arr_axpby(v[i], cb, 3, bottomBaseSize, 1 - bottomBaseSize);
		}
	}

	void cacheClear()
	{
		cache_faceIndex = -1;
		arr_set(vector_cache[0], 3 * cache_vectors_count, 0);
	}

	void cacheUpdate(size_t k)
	{
		cacheClear();

		cache_faceIndex == k;

		// Sum the bottom vertices into topCenter.
		getFace(faceVs[0], faceVs[1], faceVs[2], k);
		arr_xpy(topCenter, faceVs[0], 3);
		arr_xpy(topCenter, faceVs[1], 3);
		arr_xpy(topCenter, faceVs[2], 3);

		// Find the bottom center, the top center and an inside point.
		arr_axpy(botCenter, topCenter, 3, 1. / 3.); // base center
		scaleto(topCenter, 3, R + A); // top center (spike tip)
		double insideCoeffBottom = (useTopAsOrigin ? 0.25 : 0.75); // * bottomBaseSize / (bottomBaseSize + topBaseSize);
		arr_axpy(insidePoint, botCenter, 3, insideCoeffBottom);
		arr_axpy(insidePoint, topCenter, 3, 1 - insideCoeffBottom);

		// Copy face to bottom base
		arr_cpy(basesVs[0], faceVs[0], 9);

		// Translate the bottom face to the top and scale it.
		for(int i = 0; i < 3; i++)
			for(int dim = 0; dim < 3; dim++)
				basesVs[i + 3][dim] = topCenter[dim] + topBaseSize * (basesVs[i][dim] - botCenter[dim]);

		// Scale the bottom face
		for(int i = 0; i < 3; i++)
			arr_axpby(basesVs[i], botCenter, 3, 1 - bottomBaseSize, bottomBaseSize);

		// Find the enveloping rectangle for all vertices
		arr_cpy(rectL, basesVs[0], 3);
		arr_cpy(rectR, basesVs[0], 3);
		for(int i = 1; i < 6; i++)
			rectUnion(rectL, rectR, basesVs[i], 3);

		// Determine a reference (local origin) vertex for each face.
		// Compute the cofactors for the triple product.
		for(int f = 0; f < 5; f++)
		{
			//Select one of the faces, i.e. 3 of the six vertices.
			double *v0, *v1, *v2;
			double a[3], b[3];

			switch(f)
			{
			case 0: // Bottom
				if(!hasBottomBase)
					continue;
				v0 = basesVs[0];
				v1 = basesVs[1];
				v2 = basesVs[2];
				break;
			case 1: // top face
				if(!hasTopBase)
					continue;
				v0 = basesVs[3];
				v1 = basesVs[4];
				v2 = basesVs[5];
				break;
			default: // side face
				v0 = basesVs[f - 2 + 3 * useTopAsOrigin]; //1
				v1 = basesVs[(f - 1) % 3 + 3 * useTopAsOrigin]; //2
				v2 = basesVs[f - 2 + 3 * (1 - useTopAsOrigin)]; //4
			}

			// Find the edge vectors
			arr_cpy(a, v1, 3);
			arr_cpy(b, v2, 3);
			arr_axpy(a, v0, 3, -1);
			arr_axpy(b, v0, 3, -1);

			// Update the cache for the face origin and _cofactors
			arr_cpy(origins[f], v0, 3);
			cofactors[f][0] = a[1] * b[2] - a[2] * b[1];
			cofactors[f][1] = a[2] * b[0] - a[0] * b[2];
			cofactors[f][2] = a[0] * b[1] - a[1] * b[0];
		}

		// Compute triple product with the inside point
		for(int f = 0; f < 5; f++)
		{
			// Vector form the face origin to the inside point.
			double abc, c[3];
			arr_cpy(c, insidePoint, 3);
			arr_axpy(c, origins[f], 3, -1);

			//Triple product with the edge vectors
			abc = arr_dotprod(c, cofactors[f], 3);
			if(abc < 0)
				arr_ax(cofactors[f], 3, -1);
		}
	}

	void cacheUpdate2(size_t k)
	{
		arr_set(cache_vectors[0], 3 * cache_vectors_count, 0);
		double v[9][3]; //temp vectors
		double *bottomCenter = (double*) (v + 6); // bottom center
		double *topCenter = (double*) (v + 7); // top center (spike tip)
		double *insidePoint = (double*) (v + 8); // body (spike) center
		arr_set((double*) (v[0]), 30, 0);

		cache_faceIndex == k;

		getFace(v[0], v[1], v[2], k);

		// Sum the bottom vertices into topCenter
		for(int i = 0; i < 3; i++)
			arr_xpy(topCenter, v[i], 3);

		arr_axpy(bottomCenter, topCenter, 3, 1. / 3.); // base center
		scaleto(topCenter, 3, R + A); // top center (spike tip)
		double insideCoeffBottom = (useTopAsOrigin ? 1 : 3) * bottomBaseSize / (bottomBaseSize + topBaseSize) / 4;
		arr_axpy(insidePoint, bottomCenter, 3, insideCoeffBottom);
		arr_axpy(insidePoint, topCenter, 3, 1 - insideCoeffBottom); // an inside point

		// Translate the bottom face to the top and scale it.
		for(int i = 0; i < 3; i++)
			for(int dim = 0; dim < 3; dim++)
				v[i + 3][dim] = topCenter[dim] + topBaseSize * (v[i][dim] - bottomCenter[dim]);

		// Scale the bottom face
		for(int i = 0; i < 3; i++)
			arr_axpby(v[i], bottomCenter, 3, bottomBaseSize, 1 - bottomBaseSize);

		// Find the enveloping rectangle for all vertices
		arr_cpy(_rect[0], v[0], 3);
		arr_cpy(_rect[1], v[0], 3);
		for(int i = 1; i < 6; i++)
			rectUnion(_rect[0], _rect[1], v[i], 3);

		int face_begin = 1; // Don't check the bottom face
		int face_end = 3 + hasTopBase;

		// Determine a reference (local origin) vertex for each face.
		// Compute the cofactors for the triple product.
		for(int f = face_begin; f <= face_end; f++)
		{
			//Select one of the faces, i.e. 3 of the six vertices.
			double a[3], b[3];
			int i0, i1, i2;
			switch(f)
			{
			case 0: // Bottom
				if(!hasBottomBase)
					continue;
				i0 = 0;
				i1 = 1;
				i2 = 2;
				break;
			case 1: // side face
				i0 = 1;
				i1 = 2;
				i2 = 4;
				break;
			case 2: // side face
				i0 = 2;
				i1 = 0;
				i2 = 5;
				break;
			case 3: // side face
				i0 = 0;
				i1 = 1;
				i2 = 3;
				break;
			case 4: // top face
				if(!hasTopBase)
					continue;
				i0 = 3;
				i1 = 4;
				i2 = 5;
				break;
			}

			if(useTopAsOrigin)
			{
				i0 += 3;
				i1 += 3;
				i2 -= 3;
			}

			// Find the edge vectors
			arr_cpy(a, v[i1], 3);
			arr_cpy(b, v[i2], 3);
			arr_axpy(a, v[i0], 3, -1);
			arr_axpy(b, v[i0], 3, -1);

			// Update the cache for the face origin and _cofactors
			arr_cpy(_faceOrigins[f], v[i0], 3);
			_cofactors[f][0] = a[1] * b[2] - a[2] * b[1];
			_cofactors[f][1] = a[2] * b[0] - a[0] * b[2];
			_cofactors[f][2] = a[0] * b[1] - a[1] * b[0];
		}

		// Compute triple product with the inside point
		for(int f = face_begin; f <= face_end; f++)
		{
			// Vector form the face origin to the inside point.
			double abc, c[3];
			arr_cpy(c, insidePoint, 3);
			arr_axpy(c, _faceOrigins[f], 3, -1);

			//Triple product with the edge vectors
			abc = arr_dotprod(c, _cofactors[f], 3);
			if(abc < 0)
				arr_ax(_cofactors[f], 3, -1);
		}
	}

	double *getSpikeRectEdges(size_t k)
	{
		if(k != cache_faceIndex)
			cacheUpdate(k);

		return _rect[0];
	}

	void getSpikeEnclosingRectangle(size_t k, double** leftEdge, double** rightEdge)
	{
		if(k != cache_faceIndex)
			cacheUpdate(k);

		*leftEdge = this->rectL;
		*rightEdge = this->rectR;
	}

	bool isInSpike(size_t k, double xyz[3])
	{
		if(k != cache_faceIndex)
			cacheUpdate(k);

		double c[3]; // vector from reference vertex to xyz.
		int sign_2 = 1; // The sign of the first side. Must be 1 for the bases.
		int sign_f;

		for(int f = 0; f < 5; f++)
		{
			if((f == 0 && !hasBottomBase) || (f == 1 && !hasTopBase))
				continue;

			arr_cpy(c, xyz, 3);
			arr_axpy(c, origins[f], 3, -1);
			if(f == 0 && lenl(xyz, 3) <= R)
				return false;
			if(f == 1 && (R + A <= lenl(xyz, 3)))
				return false;

			sign_f = sign(arr_dotprod(c, cofactors[f], 3));

			if(f == 2)
			{
				if((!isHourglass && sign_f == -1) || sign_f == 0)
					return false;
				sign_2 = sign_f;
			}
			else if(sign_f != sign_2)
				return false;
		}

		return true;
	}

	void fprint_cache(FILE* fptr)
	{
		fprintf(fptr, "[ %lld , [ ", cache_faceIndex);
		for(int i = 0; i < cache_vectors_count; i++)
		{
			fprintf(fptr, "[ %e , %e , %e ], ", vector_cache[i][0], vector_cache[i][1], vector_cache[i][2]);
			switch(i)
			{
			case 2:
			case 5:
			case 8:
			case 13:
			case 18:
			case 23:
				fprintf(fptr, "\n");
			}
		}
		fprintf(fptr, "]\n");
	}

	bool isInSpikeBottomless2(size_t k, double xyz[3])
	{
		if(k != cache_faceIndex)
			cacheUpdate(k);

		double c[3]; // vector from reference vertex to xyz.
		int sign_1; // The sign of the first side.

		int f = 0; // Skip the bottom
//		if(hasBottomBase)
//		{
//			arr_cpy(c, xyz, 3);
//			arr_axpy(c, _faceOrigins[f], 3, -1);
//			if(0 >= arr_dotprod(c, _cofactors[f], 3))
//				return false;
//		}

		//side face 1
		f++;
		arr_cpy(c, xyz, 3);
		arr_axpy(c, _faceOrigins[f], 3, -1);
		sign_1 = sign(arr_dotprod(c, _cofactors[f], 3));
		if((!isHourglass && sign_1 == -1) || sign_1 == 0)
			return false;

		//side face 2
		f++;
		arr_cpy(c, xyz, 3);
		arr_axpy(c, _faceOrigins[f], 3, -1);
		if(sign_1 != sign(arr_dotprod(c, _cofactors[f], 3)))
			return false;

		//side face 3
		f++;
		arr_cpy(c, xyz, 3);
		arr_axpy(c, _faceOrigins[f], 3, -1);
		if(sign_1 != sign(arr_dotprod(c, _cofactors[f], 3)))
			return false;

		//top face
		f++;
		if(hasTopBase)
		{
			arr_cpy(c, xyz, 3);
			arr_axpy(c, _faceOrigins[f], 3, -1);
			if(0 >= arr_dotprod(c, _cofactors[f], 3))
				return false;
		}

		return true;
	}

	void writePythonFile(char* filename)
	{
		FILE* fptr = fopen(filename, "w");

		fprintf(fptr, "### nupmy must be imported as np."
				"### Sample load script:"
				"#pyscriptfile = open(datapath, \'r\')\n"
				"#pyscript = pyscriptfile.read()"
				"#pyscript.close()\n"
				"#Optionally: #pyscript = pyscript.replace(\"nan\", \"np.nan\")\n"
				"#data = eval(pyscript.read())\n"
				"### Usage of data follows:"
				"#print(data.radius)"
				"#etc.\n"
				"\n");

		fprintf(fptr, "type('', (), dict(\n");

		fprintf(fptr, "# V, E, F, R = %lld, %lld, %lld, %e\n\n", nV, nE, nF, R);

		fprintf(fptr, "radius=%e ,\n\n", R);

		fprintf(fptr, "vertices=[\n");
		for(size_t i = 0; i < nV; i++)
			fprintf(fptr, "np.array([%e,%e,%e]),\n", vx[i], vy[i], vz[i]);
		fprintf(fptr, "],\n\n");

		fprintf(fptr, "edges=(\n");
		for(size_t j = 0; j < nE; j++)
			fprintf(fptr, "(%lld, %lld),\n", ev1[j], ev2[j]);
		fprintf(fptr, "),\n\n");

		fprintf(fptr, "faces=(\n");
		for(size_t k = 0; k < nF; k++)
			fprintf(fptr, "(%lld, %lld, %lld),\n", fe1[k], fe2[k], fe3[k]);
		fprintf(fptr, ")\n\n");

		fprintf(fptr, ") #end dict\n");
		fprintf(fptr, ") #end type\n\n");
		fclose(fptr);
	}
}
;

//int testTri()
//{
//	TriSphere* ts = new TriSphere(1., 1., 3);
//	ts->writePythonFile("test_tri_data.py");
//	ts->~TriSphere();
//	TRACEF("Exitting testTri()");
//	return 0;
//}
