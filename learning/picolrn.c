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
#include <errno.h>

// hyperparameters
#define NRANDS 1024

/*
	auxiliary stuff
*/

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define SQR(x) ((x)*(x))

unsigned uRidBufSize;

/*  Raw Intensity Data
	- loads an 8-bit grey image saved in the <RID> file format
	- <RID> file contents:
		- a 32-bit signed integer w (image width)
		- a 32-bit signed integer h (image height)
		- an array,pixels[], of w*h unsigned bytes representing pixel intensities
*
* load rid file data into pixels buffer.
*
* path : nonfaces/00000-IMG_2121.JPG.rid
*/
int loadrid(uint8_t* pixels[], int* nrows, int* ncols, const char* path)
{
	FILE* file;
	int w, h;

	// open file
	file = fopen(path, "rb");

	if(!file)
	{
		printf("%s opening failure %d\n", path, errno);
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
		printf("malloc error %d\n",errno);
		fclose(file);
		return 0;
	}
	uRidBufSize += w*h*sizeof(uint8_t);//total allocation buffer size
	// read image data
	fread(*pixels, sizeof(uint8_t), w*h, file);

	// clean up
	fclose(file);

	// we're done
	printf("memory used : %u\n", uRidBufSize);
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

	/* randomly genereated 2-pixel coordination to be used to compare pixel intensity.
	 * 4 bytes: two (x,y) position pairs to be compaired, each byte is 0-255.
	 * the dim of rid is normalized to 0-1.0 mapping to 0-255 coodination.
	 * size of tcodes : (int32_t*)malloc(((1<<d)-1)*sizeof(int32_t));
	 * each element, tcodes[n], is a decion feature in a node of a tree.
	 * tcodes[] : The tree has depth d, so the full complete binary tree has (2^d -1) nodes.
	 */
	int32_t* tcodes;

	 /* The level at d+1 has 2^d nodes. lut[] contains these 2^d nodes.
	  * These nodes can be regarded as leaf nodes.
	  * these nodes are the regression value after the tree traversal from the root.
	 */
	float* lut;

} rtree;

/*
 * pixels : rid raw data
 * r,c : bounding box center at (row,col)
 * nrows, ncols : object sample's dimension of rid file
 * ldim : the column of a rid raw data, it's a row-major of 2-d pixel array
 * sr, sc : currently it's diameter of a bounding box, s, because tsr and tsc are 1.0 only.
 * tcode : randomly genereated 2-pixel pair position for later comparison to generate decision tree
 * 4 bytes (x0,y0,x1,y1), x, y is one byte which is 0-255.
*/
int bintest(int tcode, int r, int c, int sr, int sc, uint8_t pixels[],
			int nrows, int ncols, int ldim)
{
	//
	int r1, c1, r2, c2;
	int8_t* p = (int8_t*)&tcode;

	//p[0],p[1] : the first pixel
	/* an object bounding box (r,c,s) is virtually normalized to 1.0 x 1.0 or 256x256
	 * and p[0],p[1],p[2],p[3] are x,y positions normalized to -128-127, so the range is 256.
	 * now transform [-128, 127] to [-0.5,0.5] ==> p[n]/256 -> [-0.5,0.5]
	 * The real coordinate is box center (r,c) + offset transformed:
	 * a pixel coordination: <r,c> + <(p/256)*sr, (p/256)*sc>,
	 * where sc and sr is the diameter of the face bounding box.
	 * since tsr and tsc are 1.0 only, sr,sc are equal to diameter of the boundong box.
	 */
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

/* Get the regression value of the object sample, pixels, by a regression tree.
 * t : one tree {depth, tcode, lut}
 * tcode : 2-pixel coordination pair for bintest. The internal nodes of a full complete binary tree.
 * an object bounding box center at (r,c) with diameter sr=sc=s
 * pixels : rid buffer of a sample
 * nrows, ncols : rid's dimension
 * ldim : rid is row-major 2d array, ldim is the column size
 */
float get_rtree_output(rtree* t, int r, int c, int sr, int sc, uint8_t pixels[],
					   int nrows, int ncols, int ldim)
{
	int d, idx;

	//
	idx = 0;

	/*
	 * tcodes : the internal nodes of a tree. each node contains a feature to check.
	 * A feature is a two-pixel pair coordination to compare the pixel intensity.
	 *              0                  0
	 *      1              2           1
	 *   3     4       5       6       2
	 * 7  8  9  10  11  12  13  14     3
	 * the left descendent of node i is (2*i + 1)
	 * the right descendent of node i is (2*i + 2)
	 */
	for(d=0; d<t->depth; ++d)//iterate over tree depth by bintest.
	{
		if( bintest(t->tcodes[idx], r, c, sr, sc, pixels, nrows, ncols, ldim) )
			idx = 2*idx + 2;	//bintest > 0, the right descendent node of a node idx of a full binary tree
		else
			idx = 2*idx + 1;//bintest <=0, the left descendent node of a node idx of a full binary tree
		printf("(d=%d, idx=%d)\n", d, idx);
	}

	/*
	 * A full complete binary tree with depth d (from 0 - (d-1)) which has 2^d -1 nodes,
	 *
	 * lut[] : the bottom nodes of a full binary tree. lut[] contains 2^d nodes at level d. (level base is 0)
	 * lut[] : is an array for the level d nodes with float-type values.
	 * tcodes : the internal nodes of a full binary tree. Each node has the pixel pair for bintest.
	 *                         0                    d=0 , tcodes[]
	 *                 1                2           d=1 , tcodes[]
	 *             3      4        5        6       d=2 , tcodes[], tcodes[] has a depth 3 tree.
	 *           7  8   9  10   11  12   13  14     d=depth=3 ==> lut[] contains only bottom nodes with float type.
	 *           0  1   2   3    4   5    6   7     offset from 0 for lut[].
	 * the idx of the level d starts from (2^d-1).
	 * So (idx - (2^d - 1)) is index in the level d with starting offset 0.
	 */
	printf("lut[]=%f, [%d]\n", t->lut[ idx - ((1<<t->depth)-1) ], idx - ((1<<t->depth)-1));

	/* lut[] has the regression value in the bottom nodes
	 * Return the regression value for the object.
	 */
	return t->lut[ idx - ((1<<t->depth)-1) ];
}

//a binary tree of depth, the total nodes of the rtree are 1+2+4+... 2^(depth-1)= 2^depth -1= (1<<depth) -1
int allocate_rtree_data(rtree* t, int d)
{
	/* tcodes are the total internal nodes of a full complete binary tree which has depth d (from 1)
	 * totoal internal nodes are (2^d -1),
	 * each internal node has sizeof(int32_t),
	 * so the total bytes of the internal nodes are (2^d -1) bytes.
	 */
	t->tcodes = (int32_t*)malloc(((1<<d)-1)*sizeof(int32_t));

	/* The nodes of a full complete binary tree at level d (d is from 0)
	 * = 2^d nodes in level d (d is from 0).
	 * level d are leaf nodes which has float-type value.
	 */
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

/*save the regression tree
a binary tree of depth, the total nodes are 1+2+4+... 2^(depth-1)= 2^depth -1= 1<<depth -1
t->tcodes : full binary tree nodes ,
t->lut : ???

*/
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

	fread(&d, sizeof(int), 1, f);//depth

	if(!allocate_rtree_data(t, d))//allocating the buffer to  store detector to ram
		return 0;
	//the total nodes of a full binary tree having depth d. (2^d -1)
	fread(t->tcodes, sizeof(int32_t), (1<<d)-1, f);//read the data from a file
	//2^d
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
float get_split_error(int tcode, float tvals[], int rs[], int cs[],
					  int srs[], int scs[], uint8_t* pixelss[],
					  int nrowss[], int ncolss[], int ldims[],
					  double ws[], int inds[], int indsnum)
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
		if( bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]],
			pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
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
int split_training_data(int tcode, float tvals[], int rs[], int cs[], int srs[],
						int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[],
						int ldims[], double ws[], int inds[], int indsnum)
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
		while( !bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]],
			pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
		{
			if( i==j )
				break;
			else
				++i;
		}

		while( bintest(tcode, rs[inds[j]], cs[inds[j]], srs[inds[j]], scs[inds[j]],
			pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]]) )
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
		if( !bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]],
			pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
			++n0;

	//
	return n0;
}

/*
nodeidx : init from 0
d : current tree depth, init from 0
maxd : maximum tree depth to traverse
inds[] : an array of total indsum elements init as {0,1,2,3... indsnum-1}
indsnum : np + nn, total number of the training samples to be processed
*/
int grow_subtree(rtree* t, int nodeidx, int d, int maxd, float tvals[], int rs[],
				 int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[],
				 int ncolss[], int ldims[], double ws[], int inds[], int indsnum)
{
	int i, nrands;

	int32_t tcodes[2048];//randomly generated 2-pixel coordinations.
	/* tcodes[] : Every 2-pixel pair is a feature which is a decision maker.
	 * It will split the object sample to left or right subtree.
	 * All samples passing the decision maker go into two subtrees.
	 * WMSE is computed by the two subtrees classified by the feature/hypothesis.
	 * We will select the minimum WMSE by the feature/hypothesis in tcodes[];
	 * The selected feature/hypothesis is used in the internal node.
	 * The samples going to two subtrees are recursively classified till the max depth.
	 */
	float spliterrors[2048], bestspliterror;

	int n0;

	//in order to speed up the learning, so a bounded tree depth is used to stop the recursive.
	//maxd(max tree depth) is specified by python script to start the train.
	if(d == maxd)
	{
		int lutidx;
		double tvalaccum, wsum;

		//lutidx : offset from 0 for the bottom level.
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
	else if(indsnum <= 1)//indsnum=(np + nn):total number of the training samples to be processed
	{	//only one sample in the split, so terminates the split.
		//terminal node
		t->tcodes[nodeidx] = 0;

		/* nodeidx : the node index of a full complete binary tree.
		 * The left chid		->	2*nodeidx + 1
		 * The right child	->	2*nodeidx + 2
		 *
		 */
		grow_subtree(t, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, srs, scs, pixelss,
					 nrowss, ncolss, ldims, ws, inds, indsnum);//left sub tree
		grow_subtree(t, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, srs, scs, pixelss,
					 nrowss, ncolss, ldims, ws, inds, indsnum);//right sub tree

		return 1;
	}

	/* generate a position list. each has a 2-pixel pair (4 bytes) to
	 * compare the pixel intensity.
	 * This is very similar to BRIEF keypoint descriptor.
	 */
	nrands = NRANDS;	//total 1024 2-pixel pairs are generated randomly.
					//This 2-pixel pair will be an internal node in decision tree.

	/* How to calculate Weighted Mean Square Error to generate a decision tree?
	 * Find the best attribute/feature to split the training set.
	 * Attribute/feature of an image is a two-pixel comparision.
	 * So there are many possible splits of a training set by the any
	 * two-pixel pair, because there are many two-pixel pair in an image.
	 * The training samples/set is divided into two classes/groups(+1,-1)
	 * by each feature generated from generate_binary_test()
	 * the split is a "hypothesis" to face(+1) and non-face (-1),
	 * so it may be different to the ground-truth.
	 * Each feature will be a root node of binary stump :
	 * http://en.wikipedia.org/wiki/Decision_stump
	 * Each feature then splits the whole training set into two subtrees,
	 * then the Weighted Mean Sqaure Error can be calculated.
	 * Then find the minimum WSME from all the calculated WSME for each node(feature).
	 * The feature which has the minumum WMSE is the root node to split the training set.
	 * This recursive iteration from tree top to the desired depth.
	 * It's not possible to use all possible features to generate the decision tree.
	 */
	for(i=0; i<nrands; ++i)	//total 1024 pairs(features) to be used for split the training sample
		tcodes[i] = mwcrand();//get the randomly genereated 2-pixel pair position for later comparison to generate decision tree

	#pragma omp parallel for
	for(i=0; i<nrands; ++i)//compute every WMSE by the 2-pixel pairs stored in tcodes array.
		spliterrors[i] = get_split_error(tcodes[i], tvals, rs, cs, srs, scs,
										 pixelss, nrowss, ncolss, ldims, ws,
										inds, indsnum);

	//the spliterrors array contains every WMSE calculated by the binary_test pixel pair.
	bestspliterror = spliterrors[0];
	t->tcodes[nodeidx] = tcodes[0];//???

	for(i=1; i<nrands; ++i)//get the minimum WMSE as the best attribute to split the training sample to build decision tree.
		if(bestspliterror > spliterrors[i])
		{
			bestspliterror = spliterrors[i];
			t->tcodes[nodeidx] = tcodes[i];//record the best split attribute (2-pixel pair)
		}

	//split the training data into two sub-trees by the minimum WSME attribute/feature
	//which becomes an internal node of the decision tree.
	n0 = split_training_data(t->tcodes[nodeidx], tvals, rs, cs, srs, scs, pixelss,
							 nrowss, ncolss, ldims, ws, inds, indsnum);

	//now we get two subtrees after the split. This is recursive iteration to split the two sub-trees.
	grow_subtree(t, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, srs, scs, pixelss,
				 nrowss, ncolss, ldims, ws, &inds[0], n0);
	grow_subtree(t, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, srs, scs, pixelss,
				 nrowss, ncolss, ldims, ws, &inds[n0], indsnum-n0);

	//
	return 1;
}

/*
n : np + nn, total number of the training samples
*/
int grow_rtree(rtree* t, int d, float tvals[], int rs[], int cs[], int srs[],
			   int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[],
			   int ldims[], double ws[], int n)
{
	int i,ret=0;
	int* inds;

	if(!allocate_rtree_data(t, d))
		return ret;

	//inds[] : an array which has n=(np + nn) elements.
	if(inds = (int*)malloc(n*sizeof(int))){
		//inds[]={0,1,2,3,...,n-1}
		for(i=0; i<n; ++i)
			inds[i] = i;//0,1,2,3,4

		//if(!grow_subtree(t, 0, 0, d, tvals, rs, cs, srs, scs, pixelss, nrowss,
		//				ncolss, ldims, ws, inds, n))
		ret = grow_subtree(t, 0, 0, d, tvals, rs, cs, srs, scs, pixelss, nrowss,
						ncolss, ldims, ws, inds, n);
		free(inds);
		return !!ret;
	}
	return ret;
}

/*
	object samples
*/

#define MAXNUMOS 1000000

unsigned uTotalRIDFiles;

static int numos;
/* Object bounding box center : (r,c,s)
 * r : row of the center
 * c : column of the center
 * s : diameter of the center
 */
static float ors[MAXNUMOS];
static float ocs[MAXNUMOS];
static float oss[MAXNUMOS];

//rid's raw data :8 bit grey level
static uint8_t* opixelss[MAXNUMOS];
static int onrowss[MAXNUMOS];
static int oncolss[MAXNUMOS];

/*
 * loading face rid into memory
 * folder : "faces" folder containing facennn.rid and list.txt
 * faces/list.txt : the name list of the rid file
 * face0.rid
	7	r=1.864	r=51	c=51	s=68
	53.368931 52.498096 66.367462
	50.208366 50.876999 62.476379
	49.590645 52.140343 70.164223
	47.852539 51.611729 63.247417
	48.665779 51.724995 63.176521
	54.229324 48.109943 63.397022
	48.860561 49.702484 67.248596
 * face1.rid .....
 * RID file: Raw Intensity Data
	- loads an 8-bit grey image saved in the <RID> file format
	- <RID> file contents:
		- a 32-bit signed integer w (image width)
		- a 32-bit signed integer h (image height)
		- an array,pixels[], of w*h unsigned bytes representing pixel intensities
 */
int load_object_samples(const char* folder)
{
	char buffer[1024], tempstr[100];
	FILE* list;
	//
	printf("Loading object samples from '%s'\n", folder);

	//load faces/list.txt
	sprintf(buffer, "%s/%s", folder, "list.txt");
	printf("opening %s...\n", buffer);
	list = fopen(buffer, "r");

	if(!list){
		printf("%s opening failure\n", buffer);
		return 0;
	}
	//number of object sample
	numos = 0;
	//get the total rid files number.
	/*if(fgets(tempstr, 100, list) != NULL){
		sscanf(tempstr, "%u", &uTotalRIDFiles);
		printf("total %u face rid files\n", uTotalRIDFiles);
	*/
	if(fscanf(list, "%s\n", tempstr) == 1){
		sscanf(tempstr, "%u", &uTotalRIDFiles);
		printf("total %u face rid files\n", uTotalRIDFiles);
	}else{
		printf("%s reading first line error %d\n", buffer, errno);
		return 0;
	}
	// read an rid file by the list.txt, facennnn.rid, every iteration.
	//while(fscanf(list, "%s", buffer) == 1 && (numos < MAXNUMOS))
	while(fscanf(list, "%s\n", buffer) == 1 && (numos < MAXNUMOS))
	{
		char fullpath[1024];

		int nrows, ncols;
		uint8_t* opixels;
		int i, n;

		// load rid image, face0.rid, into memory opixels[]
		sprintf(fullpath, "%s/%s", folder, buffer);
		printf("loading [%s]\n", fullpath);
		/*raw intensity data :
		 * 4 bytes width : x : ncols : rid's cols
		 * 4 bytes height : y : nrows : rid's rows
		 * opixels : 8 bit grey data loaded
		*/
		if(!loadrid(&opixels, &nrows, &ncols, fullpath))
			return 0;

		// number of samples associated with this image
		// usually "7" face bounding boxes for a rid file.
		//if(fscanf(list, "%s\n", tempstr) != 1){
		if(fgets(tempstr, 100, list) == NULL){
			printf("read string error %d\n", errno);
			return 0;
		}
		printf("+++[%s]\n",tempstr);
		int temp=0;
		if(temp = sscanf(tempstr, "%d", &n) != 1){
			printf("xxx[%d] n=%d\n",temp, n);
			return 0;
		}
		/*
		if(fscanf(list, "%d", &n) != 1){
			printf("xxx n=%d\n", n);
			return 0;
		}*/

		/* One rid file buffer, opixels, will be used by
		 * "7" face bounding boxes. These bounding boxes
		 * are a little variant in the r,c,s.
		 *
		 * The total object samples are multiplied by 7 times.
		 * uTotalRIDFiles * 7
		 */
		for(i=0; i<n; ++i)
		{
			//face bounding box is specified by three coordinates:
			//(r,c,s) = (row of box center, column of box center , diameter of the box)
			float r, c, s;

			//
			if(fscanf(list, "%f %f %f", &r, &c, &s) != 3){
				printf("xxxx : r,c,s\n");
				return 0;
			}
			//object bounding box center coordinations:(r,c,s)
			ors[numos] = r;
			ocs[numos] = c;
			oss[numos] = s;

			opixelss[numos] = opixels;//8 bit raw intensity data
			onrowss[numos] = nrows;	//image rows
			oncolss[numos] = ncols;	//image cols

			/* total object samples are loaded.
			 * The total object samples are multiplied by 7 times.
			 * uTotalRIDFiles * 7
			 */
			++numos;
		}
	}
	fclose(list);
	printf("+++total faces numos=%d\n", numos);
	//check the maximum numbers of object samples supported.
	/*if(numos >= MAXNUMOS)
	{
		printf("maximum allowed number of object samples exceeded: terminating ...\n");
		return 0;
	}*/

	return 1;
}

/*
	background images
*/

#define MAXNUMBS 100000
unsigned uTotalBSRIDFiles;
static int numbs = 0;

static uint8_t* bpixelss[MAXNUMBS];
int bnrowss[MAXNUMBS];
int bncolss[MAXNUMBS];
unsigned ubgRidBufSize;

/* loading nonfaces : nonfaces/list.txt
 * 00000-IMG_2121.JPG.rid
 *
 */
int load_background_images(char* folder)
{
	FILE* list;
	char path[1024], name[1024];

	//
	printf("\nLoading background images from '%s'\n", folder);

	//
	sprintf(path, "%s/list.txt", folder);

	list = fopen(path, "r");
	printf("loading [%s]\n", path);
	if(!list)
		return 0;

	//get the total rid files number.
	if(fscanf(list, "%u\n", &uTotalBSRIDFiles) == 1){
		printf("total [%u] nonface rid files in the list...\n", uTotalBSRIDFiles);
	}else{
		printf("%s reading first line error %d\n", path, errno);
		return 0;
	}
	//
	while(fscanf(list, "%s\n", name)==1 && numbs<MAXNUMBS)
	{
		//rid file name in nonfaces : 00000-IMG_2121.JPG.rid
		sprintf(path, "%s/%s", folder, name);
		printf("[%s]\n", path);
		//
		if(loadrid(&bpixelss[numbs], &bnrowss[numbs], &bncolss[numbs], path))
			++numbs;
	}

	//
	fclose(list);
	printf("total background [%d]\n", numbs);
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

/*
 * classifying an object whether it's an object(+1) or not (-1)
 * This is a cascade topology. The negative sample is rejected first and
 * the positive will pass to the next stage to be classified again.
o  : The regression value of an object sample.
r : object sample center at row
c : object sample center at column
s : object sample center's diameter
pixels : object sample's rid data buffer
nrows, ncols : object sample's dim
ldim : object samples' column ???
*/
int classifyObject(float* o, float r, float c, float s, uint8_t pixels[],
					int nrows, int ncols, int ldim)
{
	int i, j;
	float ir, ic, isr, isc;

	//printf("[%s]: odetector.numstages=%d\n", __func__, odetector.numstages);
	//regression value of a sample image.
	*o = 0.0f;

	//0 stage. it's in init phase only!
	if(!odetector.numstages)
		return 1;

	//
	ir = (int)( r );
	ic = (int)( c );

	//currently odetector.tsr = odetector.tsc = 1.0.
	isr = (int)( odetector.tsr*s );//=s, object sample center's diameter
	isc = (int)( odetector.tsc*s );//=s, object sample center's diameter

	//
	i = 0;

	while(i < odetector.numstages)//cascade classification
	{	//get the sum of all the trees output from a training sample.
		//get_rtree_output : get the regression value of the object sample, pixels, by a regression tree.
		for(j=0; j<odetector.numtreess[i]; ++j)
			*o += get_rtree_output(&odetector.rtreearrays[i][j], ir, ic, isr,
								   isc, pixels, nrows, ncols, ldim);

		/* If the output sum of the trees in this stage is less than threshold,
		 * reject it to negative
		 */
		if(*o <= odetector.thresholds[i])
			return -1;//negative sample is rejected at the earlier stage

		//the positive samples are passed to the second stage to verify again.
		++i;
	}

	/* An object sample has passed all the stages and its output sum is
	 * greater than the threshold, so it's a positive sample.(+1)
	 */
	return 1;
}

//save detctor to a file
int saveDetector(char* path)
{
	int i, j;

	FILE* f = fopen(path, "wb");

	if(!f)
		return 0;

	//
	fwrite(&odetector.tsr, sizeof(float), 1, f);
	fwrite(&odetector.tsc, sizeof(float), 1, f);

	fwrite(&odetector.numstages, sizeof(int), 1, f);//odetector.numstages == 0 in the detector init stage.

	/* iterate over all stages
	 * if this is the init phase, the numstages is 0.
	 */
	for(i=0; i<odetector.numstages; ++i)
	{	//the ensemble of trees in the stage
		fwrite(&odetector.numtreess[i], sizeof(int), 1, f);

		//for every regression tree in the sta
		for(j=0; j<odetector.numtreess[i]; ++j)
			save_rtree_to_file(&odetector.rtreearrays[i][j], f);

		fwrite(&odetector.thresholds[i], sizeof(float), 1, f);
	}

	fclose(f);

	return 1;
}

//load a detctor from a file, d.
int loadDetector(char* path)
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
int learn_new_stage(int stageidx, int tdepth, float mintpr, float maxfpr,
					int maxnumtrees, int classs[], float rs[], float cs[],
					float ss[], uint8_t* pixelss[], int nrowss[],
					int ncolss[], float os[], int np, int nn)
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

		grow_rtree(&odetector.rtreearrays[stageidx][numtrees], tdepth, tvals,
				   irs, ics, isrs, iscs, pixelss, nrowss, ncolss, ncolss,
					ws, np+nn);

		printf("\r");
		printf("	tree %d (%f [sec]) ...", numtrees+1, getticks()-t);

		++numtrees;

		// update outputs ...
		for(i=0; i<np+nn; ++i)
		{
			float o;

			//
			o = get_rtree_output(&odetector.rtreearrays[stageidx][numtrees-1],
								 irs[i], ics[i], isrs[i], iscs[i], pixelss[i],
								nrowss[i], ncolss[i], ncolss[i]);

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
 * classifying positive(object) samples and negative(non-object) samples
 * return estimated false positive rate : FP/(TN+FP)
 *
numos : total object samples loaded. each rid sample has 7 face boxes.
each face sample has a mirrored sample, so a face sample will generate 14 face boxes, numos.

the following are init before invoking TrainedBySampleData.
#define MAXMAXNUMSAMPLES 2*MAXNUMOS	: twice the loaded samples.

float rs[MAXMAXNUMSAMPLES];	face bounding box center, row
float cs[MAXMAXNUMSAMPLES];	face bounding box center, column
float ss[MAXMAXNUMSAMPLES];	face bounding box center, diameter
int classs[MAXMAXNUMSAMPLES];	ground truth table for positive sample, +1
uint8_t* pixelss[MAXMAXNUMSAMPLES];	rid raw data, 8 bit grey level
int nrowss[MAXMAXNUMSAMPLES];	rid file's dimension
int ncolss[MAXMAXNUMSAMPLES];	rid file's dimension
float os[MAXMAXNUMSAMPLES];		regression value
int maxn;	maxnumsamples, "2*numos" numbers of object samples loaded into ram
int* np, int* nn; number of positive and negative samples
*/
float TrainedBySampleData(int classs[], float rs[], float cs[], float ss[],
						   uint8_t* pixelss[], int nrowss[], int ncolss[],
						   float os[], int maxn, int* np, int* nn)
{
	int i, n;

	int64_t nw;
	float etpr, efpr;

	int t;
	//to init the random
	#define NUMPRNGS 1024
	static int prngsinitialized = 0;
	static int64_t prngs[NUMPRNGS];

	int stop;

	//
	printf("- sampling training data (randomized) ...\n");

	t = getticks();

	//
	n = 0;
	printf("[%s]: odetector.numstages=%d\n", __func__, odetector.numstages);

	/* ===============================
	 * positive (face) object samples
	 * ===============================
	 * numos : numbers of bounding boxes loaded into ram
	 * ors, ocs and oss have been setup in load_object_samples.
	 * ors[] = object bounding box center at row array
	 * ocs[] : object bounding box center at column array
	 * oss[] : object bounding box center at diameter array
	 * opixelss : face rid raw data, 8 bit grey
	 * onrowss : rid's row
	 * oncolss : rid's column
	 * os[]: regression value of the object sample
	 */
	for(i=0; i<numos; ++i)//classifying an object i whether it's an object(+1) or not (-1).
		//os[n] : The regression value of object n after the classifier.
		if( classifyObject(&os[n], ors[i], ocs[i], oss[i], opixelss[i],
			onrowss[i], oncolss[i], oncolss[i])>0 )
		{	//This sample i is a positive sample after checking the decision tree.
			//store sample's bounding box center as object n.
			rs[n] = ors[i];
			cs[n] = ocs[i];
			ss[n] = oss[i];
			//store its rid raw file and raw file's dim
			pixelss[n] = opixelss[i];
			nrowss[n] = onrowss[i];
			ncolss[n] = oncolss[i];

			classs[n] = +1;//the sample is an object, +1

			//
			++n;
		}

	*np = n;//number of positive samples

	/*
	 * =======================================
	 * negative(non-face) non-object samples
	 * =======================================
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
	nw = 0;	//total negative
	*nn = 0;//false positive

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

			//background object
			pixels = bpixelss[idx];
			nrows = bnrowss[idx];
			ncols = bncolss[idx];

			r = mwcrand_r(&prngs[thid])%bnrowss[idx];
			c = mwcrand_r(&prngs[thid])%bncolss[idx];
			s = mwcrand_r(&prngs[thid])%( 2*MIN(MIN(r, nrows-r), MIN(c, ncols-c)) + 1 );

			if(s<24)
				continue;

			//classifying an object whether it's an object(+1) or not (-1).
			if( classifyObject(&o, r, c, s, pixels, nrows, ncols, ncols)>0 )
			{
				//we have a false positive : FPR = false_postive / (false_positive + true_nagative)
				#pragma omp critical
				{
					if(n<maxn)
					{	//store bounding box of the nonobject
						rs[n] = r;
						cs[n] = c;
						ss[n] = s;
						//store rid data and its dim of the nonobject
						pixelss[n] = pixels;
						nrowss[n] = nrows;
						ncolss[n] = ncols;

						os[n] = o;

						classs[n] = 0;//the sample is an non-object, 0

						//
						++n;//total samples of object and non-object
						++*nn;//number of false positive samples
					}
					else
						stop = 1;
				}
			}//else it's true negative.

			#pragma omp atomic
			++nw;//total negative
		}
	}

	/*
		print estimated true positive and false positive rates : false_postive / (false_positive + true_nagative)
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
sscanf(argv[5], "%f", &targetfpr);	//false positive rate : false_postive / (false_positive + true_nagative)
sscanf(argv[7], "%f", &minstagetpr);	//minimum true positive rate :
sscanf(argv[8], "%f", &maxstagefpr);	//maximum false positive rate : 0.5 is maximum because a random guess possibility is 0.5.
					//A weak classifier should be better than or equal to a random guess. the fp objects will go to the next stage
sscanf(argv[6], "%d", &tdepths);	//maximum decision tree depth. deeper gets slow response.
sscanf(argv[9], "%d", &maxnumtreesperstage);	//a tree is a weak classifier, a stage, an ensemble, is a strong classiffier with combining trees

src : path to read the the detector file "d"
dst : path to write the detector file "d"

picolrn d faces nonfaces 1 1e-6 6 0.980 0.5 1 d >> log.txt
maxnstages = maxnumstagestoappend = 1,
targetfpr=1e-6,
...
*/
int LearnByAppendStages(char* src, char* dst, int maxnumstagestoappend,
		float targetfpr, float minstagetpr, float maxstagefpr, int tdepths,
		int maxnumtreesperstage)
{
	#define MAXMAXNUMSAMPLES 2*MAXNUMOS

	static float rs[MAXMAXNUMSAMPLES];
	static float cs[MAXMAXNUMSAMPLES];
	static float ss[MAXMAXNUMSAMPLES];

	/* classs[] : positive sample or negative sample of the sampling data by a decision tree
	 * object sample is 1 and negative sample is 0
	 */
	static int classs[MAXMAXNUMSAMPLES];

	static uint8_t* pixelss[MAXMAXNUMSAMPLES];//rid
	static int nrowss[MAXMAXNUMSAMPLES];	//rid's rows
	static int ncolss[MAXMAXNUMSAMPLES];	//rid cols
	//regression value, output of sample
	static float os[MAXMAXNUMSAMPLES];

	//
	int i, maxnumsamples, maxnumstages, np, nn;

	//
	if(!loadDetector(src))
		return 0;

	/* odetector.numstages is "0" at first. It's current stage in the training phase.
	 * and then stage increases by the following learning append stage
	 */
	maxnumstages = odetector.numstages + maxnumstagestoappend;
	maxnumsamples = 2*numos;//numos : number of object samples loaded into ram

	//np : number of positive samples, nn : number of negative samples
	np = nn = 0;

	/* Only appending phase will come here, but init phase won't come here.
	 * ./picolrn d faces nonfaces '1' 1e-6 6 0.980 0.5 1 d >> log.txt
	 * i=odetector.numstages = 0 after the first init stage.
	 * maxnumstages = 0 + 1 = 1
	 * 	.....
	 */
	for(i=odetector.numstages; i<maxnumstages; ++i)
	{
		float currentfpr;//fpr= 1e-6: false positive rate, false_postive / (false_positive + true_nagative)

		printf("--------------------------------------------------------------------\n");
		printf("%d/%d\n", i+1, maxnumstages);

		/*
			sample training set
		*/

		printf("\n");
		/*
		 * classifying positive(object) samples and negative(non-object) samples.
		 * returning the estimated false positive rate : FP/(TN+FP)
		*/
		currentfpr = TrainedBySampleData(classs,		//ground truth table for positive sample, +1
														//negative, 0
										  rs, cs, ss,	//face bounding box center
										  pixelss,		//rid raw data
										  nrowss, ncolss,//rid's dim rows and cols
										os,				//regression value
										maxnumsamples,	//2*numos numbers of object samples loaded into ram
										&np, &nn);		//positive and negative samples

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

		learn_new_stage(i, tdepths, minstagetpr, maxstagefpr,
						maxnumtreesperstage, classs, rs, cs, ss,
					pixelss, nrowss, ncolss, os, np, nn);

		/*
			we have a new stage
		*/

		++odetector.numstages;

		//
		printf("\n");

		if(saveDetector(dst))
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

	//only init phase has argc==4.
	if(argc == 4)
	{
	/* create and init an object detector "d" into a file
	 * start the learning process
	 * ./picolrn 1 1 d > log.txt
	 */
		sscanf(argv[1], "%f", &odetector.tsr);
		sscanf(argv[2], "%f", &odetector.tsc);
		//
		odetector.numstages = 0;//init

		//detector d is saved to file
		if(!saveDetector(argv[3]))
			return 0;

		//
		printf("INITIALIZING: (tsr=%f, tsc=%f, numstages=%d)\n",
			   odetector.tsr, odetector.tsc, odetector.numstages);

		//
		return 0;//it only creates a "d" header file.
	}
	else if(argc == 11)
	{
	/*append stages : updating odetector.numstages
	./picolrn d faces nonfaces 1 1e-6 6 0.980 0.5 1 d >> log.txt, append 1 stage only
	./picolrn d faces nonfaces 1 1e-6 6 0.980 0.5 1 d >> log.txt, append 1 stage only
	./picolrn d faces nonfaces 1 1e-6 6 0.985 0.5 1 d >> log.txt, append 1 stage only
	./picolrn d faces nonfaces 1 1e-6 6 0.990 0.5 2 d >> log.txt, append 1 stage only
	./picolrn d faces nonfaces 1 1e-6 6 0.995 0.5 3 d >> log.txt, append 1 stage only
	./picolrn d faces nonfaces 6 1e-6 6 0.997 0.5 10 d >> log.txt, append 6 stages
	./picolrn d faces nonfaces 10 1e-6 6 0.999 0.5 20 d >> log.txt, append 10 stages
	*/
		src = argv[1];	//"d" the detector

		objspath = argv[2];	//training samples : faces images (object samples)
		nonobjimgspath = argv[3];//training samples : non-faces images (non object samples)
		//1 1e-6 6 0.980 0.5 1 d
		sscanf(argv[4], "%d", &maxnstages);
		sscanf(argv[5], "%f", &targetfpr);	//false positive rate : false_postive / (false_positive + true_nagative)
		sscanf(argv[6], "%d", &tdepths);	//maximum decision tree depth. deeper gets slow response.
		sscanf(argv[7], "%f", &minstagetpr);	//minimum true positive rate :
		sscanf(argv[8], "%f", &maxstagefpr);	//maximum false positive rate : 0.5 is maximum because a random guess possibility is 0.5.
							//A weak classifier should be better than or equal to a random guess. the fp objects will go to the next stage
		sscanf(argv[9], "%d", &maxnumtreesperstage);	//a tree is a weak classifier, a stage, an emsemble, is a strong classiffier with combining trees

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
	if(!load_object_samples(objspath))//object training samples : faces images
	{
		printf("cannot load object samples ... exiting ...\n");
		return 1;
	}
	printf("%d object samples loaded in %f [s]\n", numos, getticks()-t);

	//
	t = getticks();
	if(!load_background_images(nonobjimgspath))//non-object training samples : non-faces images
	{//load nonfaces
		printf("cannot load background images ... exiting ...\n");
		return 1;
	}
	printf("%d background images loaded in %f [s]\n", numbs, getticks()-t);

	//
	t = getticks();

	printf("LEARNING ...\n");
	/*./picolrn d faces nonfaces 1 1e-6 6 0.980 0.5 1 d >> log.txt
	 * maxnstages = 1, targetfpr=1e-6
	 *
	 */
	LearnByAppendStages(src, dst, maxnstages, targetfpr, minstagetpr,
							   maxstagefpr, tdepths, maxnumtreesperstage);
	printf("FINISHED ...\n");

	printf("elapsed time: %f [sec]\n", getticks()-t);

	printf("\n");

	//transform detector "d" to a hex array
	//python tohexarray.py d > face-detector-from-genki-dataset.ea
	//
	return 0;
}
