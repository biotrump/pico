/*
* Copyright (c) 2013, Nenad Markus
* All rights reserved.
*
* This is an implementation of the algorithm described in the following paper:
* N. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer,
* A method for object detection based on pixel intensity comparisons,
* http://arxiv.org/abs/1305.4537
*
* Redistribution and use of this program as source code or in binary form, with or without modifications, are permitted provided that the following conditions are met:
* 1. Redistributions may not be sold, nor may they be used in a commercial product or activity without prior permission from the copyright holder (contact him at nenad.markus@fer.hr).
* 2. Redistributions may not be used for military purposes.
* 3. Any published work which utilizes this program shall include the reference to the paper available at http://arxiv.org/abs/1305.4537
* 4. Redistributions must retain the above copyright notice and the reference to the algorithm on which the implementation is based on, this list of conditions and the following disclaimer.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdint.h>

// hyperparameters
#define NRANDS 1024

/*
	auxiliary stuff
*/

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define SQR(x) ((x)*(x))

/*  Raw Intensity Data
	- loads an 8-bit grey image saved in the <RID> file format
	- <RID> file contents:
		- a 32-bit signed integer w (image width)
		- a 32-bit signed integer h (image height)
		- an array,pixels[], of w*h unsigned bytes representing pixel intensities
*/

int loadrid(uint8_t* pixels[], int* nrows, int* ncols, const char* path)
{
	FILE* file;
	int w, h;

	// open file
	file = fopen(path, "rb");

	if(!file)
	{
		return 0;
	}

	// read width
	fread(&w, sizeof(int), 1, file);
	// read height
	fread(&h, sizeof(int), 1, file);

	// allocate image memory
	*nrows = h;
	*ncols = w;

	*pixels = (uint8_t*)malloc(w*h*sizeof(uint8_t));

	if(!*pixels)
	{
		fclose(file);

		return 0;
	}

	// read image data
	fread(*pixels, sizeof(uint8_t), w*h, file);

	// clean up
	fclose(file);

	// we're done
	return 1;
}

/*
	portable time function
*/

#ifdef __GNUC__
#include <time.h>
float getticks()
{
	struct timespec ts;

	if(clock_gettime(CLOCK_MONOTONIC, &ts) < 0)
	{
		printf("clock_gettime error\n");

		return -1.0f;
	}

	return ts.tv_sec + 1e-9f*ts.tv_nsec;
}
#else
#include <windows.h>
float getticks()
{
	static double freq = -1.0;
	LARGE_INTEGER lint;

	if(freq < 0.0)
	{
		if(!QueryPerformanceFrequency(&lint))
			return -1.0f;

		freq = lint.QuadPart;
	}

	if(!QueryPerformanceCounter(&lint))
		return -1.0f;

	return (float)( lint.QuadPart/freq );
}
#endif

/*
	multiply with carry PRNG
*/

uint32_t mwcrand_r(uint64_t* state)
{
	uint32_t* m;

	//
	m = (uint32_t*)state;

	// bad state?
	if(m[0] == 0)
		m[0] = 0xAAAA;

	if(m[1] == 0)
		m[1] = 0xBBBB;

	// mutate state
	m[0] = 36969 * (m[0] & 65535) + (m[0] >> 16);
	m[1] = 18000 * (m[1] & 65535) + (m[1] >> 16);

	// output
	return (m[0] << 16) + m[1];
}

uint64_t prngglobal = 0x12345678000fffffLL;

void smwcrand(uint32_t seed)
{
	prngglobal = 0x12345678000fffffLL*seed;
}

uint32_t mwcrand()
{
	return mwcrand_r(&prngglobal);
}

/*
	regression trees
*/

typedef struct
{
	int depth;
	int32_t* tcodes;

	float* lut;

} rtree;

/*
tcode : 4 bytes for (x,y) position used to compare pixel intensity
	x is 2 bytes and y is 2 bytes, too.

*/
int bintest(int tcode, int r, int c, int sr, int sc, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	//
	int r1, c1, r2, c2;
	int8_t* p = (int8_t*)&tcode;

	//p[0],p[1] : the first pixel
	r1 = (256*r + p[0]*sr)/256;
	c1 = (256*c + p[1]*sc)/256;
	//p[2],p[3] : the second pixel
	r2 = (256*r + p[2]*sr)/256;
	c2 = (256*c + p[3]*sc)/256;

	//the practical position of a pixel in the image
	r1 = MIN(MAX(0, r1), nrows-1);
	c1 = MIN(MAX(0, c1), ncols-1);

	r2 = MIN(MAX(0, r2), nrows-1);
	c2 = MIN(MAX(0, c2), ncols-1);

	//compare the intensity of the two pixels in the image, output is 0 or 1
	return pixels[r1*ldim+c1]<=pixels[r2*ldim+c2];
}

float get_rtree_output(rtree* t, int r, int c, int sr, int sc, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int d, idx;

	//
	idx = 0;

	for(d=0; d<t->depth; ++d)
		if( bintest(t->tcodes[idx], r, c, sr, sc, pixels, nrows, ncols, ldim) )
			idx = 2*idx + 2;
		else
			idx = 2*idx + 1;

	//
	return t->lut[ idx - ((1<<t->depth)-1) ];
}

int allocate_rtree_data(rtree* t, int d)
{
	//
	t->tcodes = (int32_t*)malloc(((1<<d)-1)*sizeof(int32_t));

	t->lut = (float*)malloc((1<<d)*sizeof(float));

	if(!t->tcodes || !t->lut)
	{
		free(t->tcodes);
		free(t->lut);

		t->depth = 0;

		return 0;
	}

	//
	t->depth = d;

	//
	return 1;
}

int deallocate_rtree_data(rtree* t, int d)
{
	if(t->depth)
	{
		free(t->tcodes);
		free(t->lut);

		t->tcodes = 0;
		t->lut = 0;
		t->depth = 0;
	}

	return 1;
}

int save_rtree_to_file(rtree* t, FILE* f)
{
	fwrite(&t->depth, sizeof(int), 1, f);
	fwrite(t->tcodes, sizeof(int32_t), (1<<t->depth)-1, f);
	fwrite(t->lut, sizeof(float), 1<<t->depth, f);

	return 1;
}

int load_rtree_from_file(rtree* t, FILE* f)
{
	int d;

	fread(&d, sizeof(int), 1, f);

	if(!allocate_rtree_data(t, d))
		return 0;

	fread(t->tcodes, sizeof(int32_t), (1<<d)-1, f);
	fread(t->lut, sizeof(float), 1<<d, f);

	return 1;
}

/*
split the training sample set into two subsets by the attribute/feature stored in "tcode"
which is a two-pixel pair array to compare.
Calculate the wmse of the two subset after the split by a specific feature (a two pixel intensity comparison).
refer to section 2.1 of the paper.
http://forum.biotrump.com/viewtopic.php?f=8&t=326
Please also refer to BRIEF to generate random 2 pixel pairs.

tcode : a postion of 2 pixels, 4 bytes long,(x,y), x pos is 2 bytes and y pos is 2 bytes used to compare the two pixel intensity.
tvals : ground truth table of all the training samples , -1.0f negative sample, +1.0f positive sample
inds : an array of total indsum integers init as {0,1,2,3... indsnum-1}
indsnum : np + nn, total number of the training samples to be processed

TODO : ???? the error calculation does not seem to match section 2.1.
*/
float get_split_error(int tcode, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int inds[], int indsnum)
{
	int i, j;

	double wsum, wsum0, wsum1;
	double wtvalsum0, wtvalsumsqr0, wtvalsum1, wtvalsumsqr1;

	double wmse0, wmse1;

	//init vars
	//http://forum.biotrump.com/viewtopic.php?f=8&t=326
	wsum = wsum0 = wsum1 = wtvalsum0 = wtvalsum1 = wtvalsumsqr0 = wtvalsumsqr1 = 0.0;
	/*Weighted Average:
	Sum(Wi*Vi)/Sum(Wi)
	*/
	for(i=0; i<indsnum; ++i)//for all (positive + negative )training sample images
	{
		/*TODO:
		c= bintest?1:0;
			wsum[c] += ws[inds[i]];	//sum up weight of cluster c images
			wtvalsum[c] += ws[inds[i]]*tvals[inds[i]];//sum up weighted  of ground truths in cluster c
			wtvalsumsqr[c] += ws[inds[i]]*SQR(tvals[inds[i]]);
		*/
		if( bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
		{//the image is in cluster/group 1
			wsum1 += ws[inds[i]];//sum up weight of cluster 1 images
			wtvalsum1 += ws[inds[i]]*tvals[inds[i]];//sum up weighted  of ground truths in C1
			wtvalsumsqr1 += ws[inds[i]]*SQR(tvals[inds[i]]);
		}
		else
		{//the image is in cluster/group 0
			wsum0 += ws[inds[i]];//sum up weight of cluster 0 images
			wtvalsum0 += ws[inds[i]]*tvals[inds[i]];//sum up weighted  of ground truths in C0
			wtvalsumsqr0 += ws[inds[i]]*SQR(tvals[inds[i]]);
		}

		wsum += ws[inds[i]];//Sum up all weights of the training sample images. each sample has its weight.
	}
	// ??? wsum == (wsum1+ wsum0)

	//??? This formula doesn't seem to match section 2.1 WMSE.
	wmse0 = wtvalsumsqr0 - SQR(wtvalsum0)/wsum0;//??? why divided by wsum0???
	wmse1 = wtvalsumsqr1 - SQR(wtvalsum1)/wsum1;

	//
	return (float)( (wmse0 + wmse1)/wsum );
}

/*
inds : an array of total indsum integers init as {0,1,2,3... indsnum-1}
indsnum : np + nn, total number of the training samples to be processed
*/
int split_training_data(int tcode, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int inds[], int indsnum)
{
	int stop;
	int i, j;

	int n0;

	//
	stop = 0;

	i = 0;
	j = indsnum - 1;

	while(!stop)
	{
		//
		while( !bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
		{
			if( i==j )
				break;
			else
				++i;
		}

		while( bintest(tcode, rs[inds[j]], cs[inds[j]], srs[inds[j]], scs[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]]) )
		{
			if( i==j )
				break;
			else
				--j;
		}

		//
		if( i==j )
			stop = 1;
		else
		{
			// swap
			inds[i] = inds[i] ^ inds[j];
			inds[j] = inds[i] ^ inds[j];
			inds[i] = inds[i] ^ inds[j];
		}
	}

	//
	n0 = 0;

	for(i=0; i<indsnum; ++i)
		if( !bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
			++n0;

	//
	return n0;
}

/*
nodeidx : init from 0
d : current tree depth, init from 0
maxd : maximum tree depth to traverse
inds : an array of total indsum integers init as {0,1,2,3... indsnum-1}
indsnum : np + nn, total number of the training samples to be processed
*/
int grow_subtree(rtree* t, int nodeidx, int d, int maxd, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int inds[], int indsnum)
{
	int i, nrands;

	int32_t tcodes[2048];
	float spliterrors[2048], bestspliterror;

	int n0;

	//in order to speed up the learning, so a bounded tree depth is used to stop the recursive.
	//maxd(max tree depth) is specified by python script to start the train.
	if(d == maxd)
	{
		int lutidx;
		double tvalaccum, wsum;

		//
		lutidx = nodeidx - ((1<<maxd)-1);

		// compute output: a simple average
		tvalaccum = 0.0;
		wsum = 0.0;

		for(i=0; i<indsnum; ++i)
		{
			tvalaccum += ws[inds[i]]*tvals[inds[i]];//ground truth value * weight
			wsum += ws[inds[i]];
		}

		if(wsum == 0.0)
			t->lut[lutidx] = 0.0f;
		else
			t->lut[lutidx] = (float)( tvalaccum/wsum );

		//
		return 1;
	}
	else if(indsnum <= 1)
	{	//only one sample in the split is left, so terminates the split.
		//terminal node
		t->tcodes[nodeidx] = 0;

		//????
		grow_subtree(t, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);
		grow_subtree(t, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);

		return 1;
	}

	// generate a position list. each has a 2-pixel pair (4 bytes) to compare the pixel intensity. This is very similar to BRIEF keypoint descriptor.
	nrands = NRANDS;	//total 1024 2-pixel pairs are generated randomly. This 2-pixel pair will be an internal node in decision tree.

	/*How to calculate Weighted Mean Square Error to generate a decision tree?
	Find the best attribute/feature to split the training set. Attribute/feature of an image is a two-pixel comparision.
	So there are many possible splits of a training set by the any two-pixel pair, because there are many two-pixel pair in an image.
	The training samples/set is divided into two classes/groups(+1,-1) by each feature generated from generate_binary_test()
	the split is a "hypothesis" to face(+1) and non-face (-1), so it may be different to the ground-truth.
	Each feature will be a root node of binary stump : http://en.wikipedia.org/wiki/Decision_stump
	Each feature then splits the whole training set into two subtrees, then the Weighted Mean Sqaure Error can be calculated.
	Then find the minimum WSME from all the calculated WSME for each node(feature).
	The feature which has the minumum WMSE is the root node to split the training set.
	This recursive iteration from tree top to the desired depth. It's not possible to use all possible features to generate the decision tree.
	*/
	for(i=0; i<nrands; ++i)	//total 1024 pairs(features) to be used for split the training sample
		tcodes[i] = mwcrand();//get the randomly genereated 2-pixel pair position for later comparison to generate decision tree

	#pragma omp parallel for
	for(i=0; i<nrands; ++i)//compute every WMSE by the 2-pixel pairs stored in tcodes array.
		spliterrors[i] = get_split_error(tcodes[i], tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);

	//the spliterrors array contains every WMSE calculated by the binary_test pixel pair.
	bestspliterror = spliterrors[0];
	t->tcodes[nodeidx] = tcodes[0];

	for(i=1; i<nrands; ++i)//get the minimum WMSE as the best attribute to split the training sample to build decision tree.
		if(bestspliterror > spliterrors[i])
		{
			bestspliterror = spliterrors[i];
			t->tcodes[nodeidx] = tcodes[i];//record the best split attribute (2-pixel pair)
		}

	//split the training data into two sub-trees by the minimum WSME attribute/feature
	//which becomes an internal node of the decision tree.
	n0 = split_training_data(t->tcodes[nodeidx], tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);

	//now we get two subtrees after the split. This is recursive iteration to split the two sub-trees.
	grow_subtree(t, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, &inds[0], n0);
	grow_subtree(t, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, &inds[n0], indsnum-n0);

	//
	return 1;
}

/*

n : np + nn, total number of the training samples
*/
int grow_rtree(rtree* t, int d, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int n)
{
	int i;
	int* inds;

	if(!allocate_rtree_data(t, d))
		return 0;

	//
	inds = (int*)malloc(n*sizeof(int));

	for(i=0; i<n; ++i)
		inds[i] = i;//0,1,2,3,4

	/*
	ret = grow_subtree();
	free(inds);
	return !!ret;
	*/
	if(!grow_subtree(t, 0, 0, d, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, n))
	{
		free(inds);
		return 0;
	}
	else
	{
		free(inds);
		return 1;
	}
}

/*
	object samples
*/

#define MAXNUMOS 1000000

static int numos;

static float ors[MAXNUMOS];
static float ocs[MAXNUMOS];
static float oss[MAXNUMOS];

static uint8_t* opixelss[MAXNUMOS];
static int onrowss[MAXNUMOS];
static int oncolss[MAXNUMOS];

int load_object_samples(const char* folder)
{
	char buffer[1024];
	FILE* list;

	//
	printf("Loading object samples from '%s'\n", folder);

	//
	sprintf(buffer, "%s/%s", folder, "list.txt");

	list = fopen(buffer, "r");

	if(!list)
		return 0;

	//number of object sample
	numos = 0;

	while(fscanf(list, "%s", buffer) == 1) // read an image file name
	{
		char fullpath[1024];

		int nrows, ncols;
		uint8_t* opixels;
		int i, n;

		//
		if(numos >= MAXNUMOS)
		{
			printf("maximum allowed number of object samples exceeded: terminating ...\n");

			return 0;
		}

		// load rid image into memory, opixels[]
		sprintf(fullpath, "%s/%s", folder, buffer);
		/*raw intensity data :rid
		nrows : rid's rows
		ncols : rid's cols
		opixels : 8 bit grey data loaded
		*/
		if(!loadrid(&opixels, &nrows, &ncols, fullpath))
			return 0;

		// number of samples associated with this image
		if(fscanf(list, "%d", &n) != 1)
			return 0;

		// get samples
		for(i=0; i<n; ++i)
		{
			//object sample is specified by three coordinates (row, column and size;
			//all in pixels)
			float r, c, s;

			//
			if(fscanf(list, "%f %f %f", &r, &c, &s) != 3)
				return 0;

			//object 
			ors[numos] = r;
			ocs[numos] = c;
			oss[numos] = s;

			opixelss[numos] = opixels;//8 bit raw intensity data
			onrowss[numos] = nrows;	//image rows
			oncolss[numos] = ncols;	//image cols

			//total object samples are loaded.
			++numos;
		}
	}

	fclose(list);

	return 1;
}

/*
	background images
*/

#define MAXNUMBS 100000

static int numbs = 0;

static uint8_t* bpixelss[MAXNUMBS];
int bnrowss[MAXNUMBS];
int bncolss[MAXNUMBS];

//loading nonfaces
int load_background_images(char* folder)
{
	FILE* list;
	char path[1024], name[1024];

	//
	printf("Loading background images from '%s'\n", folder);

	//
	sprintf(path, "%s/list.txt", folder);

	list = fopen(path, "r");

	if(!list)
		return 0;

	//
	while(fscanf(list, "%s", name)==1 && numbs<MAXNUMBS)
	{
		//
		sprintf(path, "%s/%s", folder, name);

		//
		if(loadrid(&bpixelss[numbs], &bnrowss[numbs], &bncolss[numbs], path))
			++numbs;
	}

	//
	fclose(list);

	//
	return 1;
}

/*

*/

struct
{
	float tsr, tsc;

	int numstages;

	float thresholds[1024];
	int numtreess[1024];

	rtree rtreearrays[32][256];

} odetector;

int classify_region(float* o, float r, float c, float s, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int i, j;
	float ir, ic, isr, isc;

	//
	*o = 0.0f;

	//
	if(!odetector.numstages)
		return 1;

	//
	ir = (int)( r );
	ic = (int)( c );

	isr = (int)( odetector.tsr*s );
	isc = (int)( odetector.tsc*s );

	//
	i = 0;

	while(i < odetector.numstages)
	{
		//
		for(j=0; j<odetector.numtreess[i]; ++j)
			*o += get_rtree_output(&odetector.rtreearrays[i][j], ir, ic, isr, isc, pixels, nrows, ncols, ldim);

		//
		if(*o <= odetector.thresholds[i])
			return -1;

		//
		++i;
	}

	//
	return 1;
}

int save_to_file(char* path)
{
	int i, j;

	FILE* f = fopen(path, "wb");

	if(!f)
		return 0;

	//
	fwrite(&odetector.tsr, sizeof(float), 1, f);
	fwrite(&odetector.tsc, sizeof(float), 1, f);

	fwrite(&odetector.numstages, sizeof(int), 1, f);//odetector.numstages == 0 in the detector init stage.

	//
	for(i=0; i<odetector.numstages; ++i)
	{
		fwrite(&odetector.numtreess[i], sizeof(int), 1, f);

		for(j=0; j<odetector.numtreess[i]; ++j)
			save_rtree_to_file(&odetector.rtreearrays[i][j], f);

		fwrite(&odetector.thresholds[i], sizeof(float), 1, f);
	}

	fclose(f);

	return 1;
}

int load_from_file(char* path)
{
	int i, j;

	FILE* f = fopen(path, "rb");

	if(!f)
		return 0;

	//
	fread(&odetector.tsr, sizeof(float), 1, f);
	fread(&odetector.tsc, sizeof(float), 1, f);

	fread(&odetector.numstages, sizeof(int), 1, f);

	//
	for(i=0; i<odetector.numstages; ++i)
	{
		fread(&odetector.numtreess[i], sizeof(int), 1, f);

		for(j=0; j<odetector.numtreess[i]; ++j)
			load_rtree_from_file(&odetector.rtreearrays[i][j], f);

		fread(&odetector.thresholds[i], sizeof(float), 1, f);
	}

	fclose(f);

	return 1;
}

/*
classs : ground truth table for all training samples, +1(positive) and 0(negative)
np : number of positive samples
nn : number of negative samples
static float rs[MAXMAXNUMSAMPLES];
static float cs[MAXMAXNUMSAMPLES];
static float ss[MAXMAXNUMSAMPLES];
*/
int learn_new_stage(int stageidx, int tdepth, float mintpr, float maxfpr, int maxnumtrees, int classs[], float rs[], float cs[], float ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], float os[], int np, int nn)
{
	int i;

	float* tvals;

	int* irs;
	int* ics;
	int* isrs;
	int* iscs;

	double* ws;
	double wsum;

	float threshold, tpr, fpr;

	int numtrees;

	//tvals : the ground truth table (floating type) for positive and negative samples
	tvals = (float*)malloc((np+nn)*sizeof(float));

	//TODO, tvals[i] = classs[i]?+1.0f:-1.0f;
	for(i=0; i<np+nn; ++i)//iterating all training samples
		if(classs[i])//+1, postivie
			tvals[i] = +1.0f; // object
		else//0, negative
			tvals[i] = -1.0f; // non-object

	//
	irs = (int*)malloc((np+nn)*sizeof(int));
	ics = (int*)malloc((np+nn)*sizeof(int));

	isrs = (int*)malloc((np+nn)*sizeof(int));
	iscs = (int*)malloc((np+nn)*sizeof(int));

	for(i=0; i<np+nn; ++i)
	{
		irs[i] = (int)( rs[i] );
		ics[i] = (int)( cs[i] );

		isrs[i] = (int)( odetector.tsr*ss[i] );
		iscs[i] = (int)( odetector.tsc*ss[i] );
	}

	//
	ws = (double*)malloc((np+nn)*sizeof(double));

	//
	numtrees = 0;

	fpr = 1.0f;

	while(numtrees<maxnumtrees && fpr>maxfpr)
	{
		float t;
		int numtps, numfps;

		// compute weights ...
		wsum = 0.0;

		for(i=0; i<np+nn; ++i)
		{
			if(classs[i])
				ws[i] = exp(-1.0*os[i])/np;
			else
				ws[i] = exp(+1.0*os[i])/nn;

			wsum += ws[i];
		}

		for(i=0; i<np+nn; ++i)
			ws[i] /= wsum;//normalize the sum of ws to 1

		// grow a tree ...
		t = getticks();

		grow_rtree(&odetector.rtreearrays[stageidx][numtrees], tdepth, tvals, irs, ics, isrs, iscs, pixelss, nrowss, ncolss, ncolss, ws, np+nn);

		printf("\r");
		printf("	tree %d (%f [sec]) ...", numtrees+1, getticks()-t);

		++numtrees;

		// update outputs ...
		for(i=0; i<np+nn; ++i)
		{
			float o;

			//
			o = get_rtree_output(&odetector.rtreearrays[stageidx][numtrees-1], irs[i], ics[i], isrs[i], iscs[i], pixelss[i], nrowss[i], ncolss[i], ncolss[i]);

			//
			os[i] += o;
		}

		// get threshold ...
		threshold = 5.0f;

		do
		{
			//
			threshold -= 0.005f;

			numtps = 0;
			numfps = 0;

			//
			for(i=0; i<np+nn; ++i)
			{
				if( classs[i] && os[i]>threshold)
					++numtps;
				if(!classs[i] && os[i]>threshold)
					++numfps;
			}

			//
			tpr = numtps/(float)np;
			fpr = numfps/(float)nn;
		}
		while(tpr<mintpr);

		printf(" tpr=%f, fpr=%f\t", tpr, fpr);
		fflush(stdout);
	}

	//
	odetector.thresholds[stageidx] = threshold;

	printf("\n");
	printf("	threshold set to %f\n", threshold);

	//
	odetector.numtreess[stageidx] = numtrees;

	//
	free(tvals);

	free(irs);
	free(ics);
	free(isrs);
	free(iscs);

	free(ws);

	//
	return 1;
}

/*
numos : total object samples loaded
*/
float sample_training_data(int classs[], float rs[], float cs[], float ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], float os[], int maxn, int* np, int* nn)
{
	int i, n;

	int64_t nw;
	float etpr, efpr;

	int t;

	#define NUMPRNGS 1024
	static int prngsinitialized = 0;
	static int64_t prngs[NUMPRNGS];

	int stop;

	//
	printf("- sampling training data (randomized) ...\n");

	t = getticks();

	//
	n = 0;

	/*
		positive (face) object samples
	*/

	for(i=0; i<numos; ++i)
		if( classify_region(&os[n], ors[i], ocs[i], oss[i], opixelss[i], onrowss[i], oncolss[i], oncolss[i])>0 )
		{
			//
			rs[n] = ors[i];
			cs[n] = ocs[i];
			ss[n] = oss[i];

			pixelss[n] = opixelss[i];
			nrowss[n] = onrowss[i];
			ncolss[n] = oncolss[i];

			classs[n] = +1;//positive sample, +1

			//
			++n;
		}

	*np = n;//number of positive samples

	/*
		negative(non-face) non-object samples
	*/

	if(!prngsinitialized)
	{
		// initialize a PRNG for each thread
		for(i=0; i<NUMPRNGS; ++i)
			prngs[i] = 0xFFFF*mwcrand() + 0xFFFF1234FFFF0001LL*mwcrand();

		//
		prngsinitialized = 1;
	}

	//
	nw = 0;
	*nn = 0;

	stop = 0;

	#pragma omp parallel
	{
		int thid;

		//
		thid = omp_get_thread_num();

		while(!stop)
		{
			float o;
			int idx, s, r, c, nrows, ncols;
			uint8_t* pixels;

			//
			idx = mwcrand_r(&prngs[thid])%numbs;

			//
			pixels = bpixelss[idx];
			nrows = bnrowss[idx];
			ncols = bncolss[idx];

			r = mwcrand_r(&prngs[thid])%bnrowss[idx];
			c = mwcrand_r(&prngs[thid])%bncolss[idx];
			s = mwcrand_r(&prngs[thid])%( 2*MIN(MIN(r, nrows-r), MIN(c, ncols-c)) + 1 );

			if(s<24)
				continue;

			//
			if( classify_region(&o, r, c, s, pixels, nrows, ncols, ncols)>0 )
			{
				//we have a false positive ...
				#pragma omp critical
				{
					if(n<maxn)
					{
						rs[n] = r;
						cs[n] = c;
						ss[n] = s;

						pixelss[n] = pixels;
						nrowss[n] = nrows;
						ncolss[n] = ncols;

						os[n] = o;

						classs[n] = 0;//negative sample

						//
						++n;
						++*nn;////number of negative samples
					}
					else
						stop = 1;
				}
			}

			#pragma omp atomic
			++nw;
		}
	}

	/*
		print estimated true positive and false positive rates
	*/

	etpr = *np/(float)numos;
	efpr = (float)( *nn/(double)nw );

	printf("	tpr (sampling): %.8f\n", etpr);
	printf("	fpr (sampling): %.8f (%d/%lld)\n", efpr, *nn, (long long int)nw);
	printf("	elapsed time: %f  [sec]\n", getticks()-t);

	/*

	*/

	return efpr;
}

/*
maxnumstagestoappend : 		sscanf(argv[4], "%d", &maxnstages);
sscanf(argv[5], "%f", &targetfpr);	//false positive rate : error postive / total detected objects
sscanf(argv[7], "%f", &minstagetpr);	//minimum true positive rate : 
sscanf(argv[8], "%f", &maxstagefpr);	//maximum false positive rate : 0.5 is maximum because a random guess possibility is 0.5.
					//A weak classifier should be better than or equal to a random guess. the fp objects will go to the next stage
sscanf(argv[6], "%d", &tdepths);	//maximum decision tree depth. deeper gets slow response.
sscanf(argv[9], "%d", &maxnumtreesperstage);	//a tree is a weak classifier, a stage, en emsemble, is a strong classiffier with combining trees
*/
int append_stages_to_odetector(char* src, char* dst, int maxnumstagestoappend, float targetfpr, float minstagetpr, float maxstagefpr, int tdepths, int maxnumtreesperstage)
{
	#define MAXMAXNUMSAMPLES 2*MAXNUMOS

	static float rs[MAXMAXNUMSAMPLES];
	static float cs[MAXMAXNUMSAMPLES];
	static float ss[MAXMAXNUMSAMPLES];

	static int classs[MAXMAXNUMSAMPLES];

	static uint8_t* pixelss[MAXMAXNUMSAMPLES];
	static int nrowss[MAXMAXNUMSAMPLES];
	static int ncolss[MAXMAXNUMSAMPLES];

	static float os[MAXMAXNUMSAMPLES];

	//
	int i, maxnumsamples, maxnumstages, np, nn;

	//
	if(!load_from_file(src))
		return 0;

	//
	maxnumstages = odetector.numstages + maxnumstagestoappend;
	maxnumsamples = 2*numos;

	//
	np = nn = 0;

	for(i=odetector.numstages; i<maxnumstages; ++i)
	{
		float currentfpr;

		printf("--------------------------------------------------------------------\n");
		printf("%d/%d\n", i+1, maxnumstages);

		/*
			sample training set
		*/

		printf("\n");

		currentfpr = sample_training_data(classs, rs, cs, ss, pixelss, nrowss, ncolss, os, maxnumsamples, &np, &nn);

		if(currentfpr <= targetfpr)
		{
			printf("- target FPR achieved ... terminating learning process ...\n");

			break;
		}

		/*
			learn decision trees in the current stage
		*/

		printf("\n");
		printf("- learning stage ...\n");
		printf("	npositives: %d, nnegatives: %d\n", np, nn);

		learn_new_stage(i, tdepths, minstagetpr, maxstagefpr, maxnumtreesperstage, classs, rs, cs, ss, pixelss, nrowss, ncolss, os, np, nn);

		/*
			we have a new stage
		*/

		++odetector.numstages;

		//
		printf("\n");

		if(save_to_file(dst))
			printf("- saving partial results to '%s' ...\n", dst);
		else
			printf("- saving results to '%s' has failed ...\n", dst);

		printf("\n");
	}

	printf("--------------------------------------------------------------------\n");

	//
	return 1;
}

/*

*/

const char* howto()
{
	return
		"TODO\n"
		/*
		"Welcome to the tool for learning object detectors based on pixel intensity comparisons organized in decision trees. The implementation closely follows the following technical report:\n"
		"\n"
		"\tN. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer.\n"
		"\tA Method for Object Detection Based on Pixel Intensity Comparisons.\n"
		"\thttp://arxiv.org/abs/1305.4537\n"
		"\n"
		"RID image container ...\n"
		"\n"
		"\n"
		"Copyright (c) 2013, Nenad Markus\n"
		"All rights reserved.\n"
		*/
	;
}

int main(int argc, char* argv[])
{
	float t;

	char* src;
	char* dst;
	char* objspath;
	char* nonobjimgspath;

	int maxnstages;
	float targetfpr;
	float minstagetpr, maxstagefpr;
	int tdepths, maxnumtreesperstage;

	//
	//create and init an object detector "d" into a file
	//start the learning process
	//./picolrn 1 1 d > log.txt
	if(argc == 4)
	{
		sscanf(argv[1], "%f", &odetector.tsr);
		sscanf(argv[2], "%f", &odetector.tsc);
		//
		odetector.numstages = 0;//init

		//detector d is saved to file
		if(!save_to_file(argv[3]))
			return 0;

		//
		printf("INITIALIZING: (%f, %f)\n", odetector.tsr, odetector.tsc);

		//
		return 0;//it only creates a "d" header file.
	}
	else if(argc == 11)
	{
	/*append stages : updating odetector.numstages
	./picolrn d faces nonfaces 1 1e-6 6 0.980 0.5 1 d >> log.txt
	./picolrn d faces nonfaces 1 1e-6 6 0.980 0.5 1 d >> log.txt
	./picolrn d faces nonfaces 1 1e-6 6 0.985 0.5 1 d >> log.txt
	./picolrn d faces nonfaces 1 1e-6 6 0.990 0.5 2 d >> log.txt
	./picolrn d faces nonfaces 1 1e-6 6 0.995 0.5 3 d >> log.txt
	./picolrn d faces nonfaces 6 1e-6 6 0.997 0.5 10 d >> log.txt
	./picolrn d faces nonfaces 10 1e-6 6 0.999 0.5 20 d >> log.txt
	*/
		src = argv[1];	//"d" the detector

		objspath = argv[2];	//training samples : faces images (object samples)
		nonobjimgspath = argv[3];//training samples : non-faces images (non object samples)
		//1 1e-6 6 0.980 0.5 1 d
		sscanf(argv[4], "%d", &maxnstages);
		sscanf(argv[5], "%f", &targetfpr);	//false positive rate : error postive / total detected objects
		sscanf(argv[6], "%d", &tdepths);	//maximum decision tree depth. deeper gets slow response.
		sscanf(argv[7], "%f", &minstagetpr);	//minimum true positive rate : 
		sscanf(argv[8], "%f", &maxstagefpr);	//maximum false positive rate : 0.5 is maximum because a random guess possibility is 0.5.
							//A weak classifier should be better than or equal to a random guess. the fp objects will go to the next stage
		sscanf(argv[9], "%d", &maxnumtreesperstage);	//a tree is a weak classifier, a stage, en emsemble, is a strong classiffier with combining trees

		dst = argv[10];	//detector d in the file
	}
	else
	{
		printf("%s", howto());
		return 0;
	}

	// initialize PRNG
	smwcrand(time(0));

	//only valid in append stages , ie. argc==11
	t = getticks();
	if(!load_object_samples(objspath))//training samples : faces images
	{
		printf("cannot load object samples ... exiting ...\n");
		return 1;
	}
	printf("%d object samples loaded in %f [s]\n", numos, getticks()-t);

	//
	t = getticks();
	if(!load_background_images(nonobjimgspath))//training samples : non-faces images
	{//load nonfaces
		printf("cannot load background images ... exiting ...\n");
		return 1;
	}
	printf("%d background images loaded in %f [s]\n", numbs, getticks()-t);

	//
	t = getticks();
	printf("LEARNING ...\n");
	append_stages_to_odetector(src, dst, maxnstages, targetfpr, minstagetpr, maxstagefpr, tdepths, maxnumtreesperstage);
	printf("FINISHED ...\n");

	printf("elapsed time: %f [sec]\n", getticks()-t);

	printf("\n");

	//transform detector "d" to a hex array
	//python tohexarray.py d > face-detector-from-genki-dataset.ea
	//
	return 0;
}
