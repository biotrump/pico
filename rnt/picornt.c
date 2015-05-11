/*
 *	Copyright (c) 2013, Nenad Markus
 *	All rights reserved.
 *
 *	This is an implementation of the algorithm described in the following paper:
 *		N. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer,
 *		Object Detection with Pixel Intensity Comparisons Organized in Decision Trees,
 *		http://arxiv.org/abs/1305.4537
 *
 *	Redistribution and use of this program as source code or in binary form, with or without modifications, are permitted provided that the following conditions are met:
 *		1. Redistributions may not be sold, nor may they be used in a commercial product or activity without prior permission from the copyright holder (contact him at nenad.markus@fer.hr).
 *		2. Redistributions may not be used for military purposes.
 *		3. Any published work which utilizes this program shall include the reference to the paper available at http://arxiv.org/abs/1305.4537
 *		4. Redistributions must retain the above copyright notice and the reference to the algorithm on which the implementation is based on, this list of conditions and the following disclaimer.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))

/*

*/
#include <stdio.h>
#include <stdlib.h>

#define	MIN_SIZE (128)
#define	MAX_SIZE (1024)

typedef unsigned char uint8_t;
typedef short int16_t;

int find_objects
		(
			float rs[], float cs[], float ss[], float qs[], int maxndetections,
			int (*run_detection_cascade)(float*, int, int, int, void*, int, int, int),
			void* pixels, int nrows, int ncols, int ldim,
			float scalefactor, float stridefactor, float minsize, float maxsize
		)
{
	float s;
	int ndetections;

	//
	ndetections = 0;
	s = minsize;

	while(s<=maxsize)
	{
		float r, c, dr, dc;

		//
		dr = dc = MAX(stridefactor*s, 1.0f);

		//
		for(r=s/2+1; r<=nrows-s/2-1; r+=dr)
			for(c=s/2+1; c<=ncols-s/2-1; c+=dc)
			{
				float q;
				int t;

				if(run_detection_cascade(&q, r, c, s, pixels, nrows, ncols, ldim) == 1)
				{
					if(ndetections < maxndetections)
					{
						qs[ndetections] = q;
						rs[ndetections] = r;
						cs[ndetections] = c;
						ss[ndetections] = s;

						//
						++ndetections;
					}
				}
			}

		//
		s = scalefactor*s;
	}

	//
	return ndetections;
}

/*

*/

float get_overlap(float r1, float c1, float s1, float r2, float c2, float s2)
{
	float overr, overc;

	//
	overr = MAX(0, MIN(r1+s1/2, r2+s2/2) - MAX(r1-s1/2, r2-s2/2));
	overc = MAX(0, MIN(c1+s1/2, c2+s2/2) - MAX(c1-s1/2, c2-s2/2));

	//
	return overr*overc/(s1*s1+s2*s2-overr*overc);
}

void ccdfs(int a[], int i, float rs[], float cs[], float ss[], int n)
{
	int j;

	//
	for(j=0; j<n; ++j)
		if(a[j]==0 && get_overlap(rs[i], cs[i], ss[i], rs[j], cs[j], ss[j])>0.3f)
		{
			//
			a[j] = a[i];

			//
			ccdfs(a, j, rs, cs, ss, n);
		}
}

int find_connected_components(int a[], float rs[], float cs[], float ss[], int n)
{
	int i, ncc, cc;

	//
	if(!n)
		return 0;

	//
	for(i=0; i<n; ++i)
		a[i] = 0;

	//
	ncc = 0;
	cc = 1;

	for(i=0; i<n; ++i)
		if(a[i] == 0)
		{
			//
			a[i] = cc;

			//
			ccdfs(a, i, rs, cs, ss, n);

			//
			++ncc;
			++cc;
		}

	//
	return ncc;
}

int cluster_detections(float rs[], float cs[], float ss[], float qs[], int n)
{
	int idx, ncc, cc;
	int a[4096];

	//
	ncc = find_connected_components(a, rs, cs, ss, n);

	if(!ncc)
		return 0;

	//
	idx = 0;

	for(cc=1; cc<=ncc; ++cc)
	{
		int i, k;

		float sumqs=0.0f, sumrs=0.0f, sumcs=0.0f, sumss=0.0f;

		//
		k = 0;

		for(i=0; i<n; ++i)
			if(a[i] == cc)
			{
				sumqs += qs[i];
				sumrs += rs[i];
				sumcs += cs[i];
				sumss += ss[i];

				++k;
			}

		//
		qs[idx] = sumqs; // accumulated confidence measure

		//
		rs[idx] = sumrs/k;
		cs[idx] = sumcs/k;
		ss[idx] = sumss/k;

		//
		++idx;
	}

	//
	return idx;
}

#include "sample/facefinder.c"
/*
 * input : grayscale image
 */
int pico_facedetection(void* frame, int width, int height,
	int maxdetect, float *frs, float *fcs, float *fss)
{
	int i, j;
//	float t;

	uint8_t* pixels;
	int nrows, ncols, ldim;

	#define MAXNDETECTIONS 2048
	int ndetections=0;
	float qs[MAXNDETECTIONS], rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS];

	void* gray = NULL;
	void* pyr[5] = {NULL, NULL, NULL, NULL, NULL};
	printf("%s: %d %d %d\n", __func__, width, height, maxdetect);
	/*
		IMPORTANT:
			* these parameters are highly specific for each detection cascade
			  (determine them experimentally)
	*/

	// * this function should be generated with picogen from a detection cascade output by picolrn
	int (*run_detection_cascade)(float*, int, int, int, void*, int, int, int)
		= run_facefinder;

	// * detection quality threshold (must be >= 0.0f)
	// * you can vary the TPR and FPR with this value
	// * if you're experiencing too many false positives, try a larger number here (for example, 7.5f)
	float qthreshold = 5.0f;

	// * how much to rescale the window during the multiscale detection process
	// * increasing this value leads to lower number of detections and higher processing speed
	// * for example, set to 1.2f if you're using pico on a mobile device
	float scalefactor = 1.1f;

	// * how much to move the window between neighboring detections
	// * increasing this value leads to lower number of detections and higher processing speed
	// * for example, set to 0.05f if you want really high recall
	float stridefactor = 0.1f;

	// * coarse image pyramid support
	// * can improve noise and aliasing problems in some applications
	// * set to 1 if pico fails to detect large objects
	int usepyr = 0;

	// perform detection with the pico library
//	t = getticks();

	if( (gray = malloc(width*height)) == NULL ){
		printf("!!!!gray malloc failure\n");
		return 0;
	}else{
		memcpy(gray, frame, width*height);
	}
	if(usepyr)
	{
#if 0
		int nd;
		//
		pyr[0] = gray;
		//pyr[1] = cvCreateImage(cvSize(frame->width/2, frame->height/2), frame->depth, 1);
		pyr[1] = malloc((width/2)*(height/2));
		//pyr[2] = cvCreateImage(cvSize(frame->width/4, frame->height/4), frame->depth, 1);
		pyr[2] = malloc((width/4)*(height/4));
		//pyr[3] = cvCreateImage(cvSize(frame->width/8, frame->height/8), frame->depth, 1);
		pyr[3] = malloc((width/8)*(height/8));
		//pyr[4] = cvCreateImage(cvSize(frame->width/16, frame->height/16), frame->depth, 1);
		pyr[4] = malloc((width/16)*(height/16));

		pixels = (uint8_t*)pyr[0];
		nrows = height;
		ncols = width;
		ldim = width;//widthStep;TODO TODO .... if this is the same width for row-major image, 8 bit gray pixel

		ndetections = find_objects(rs, cs, ss, qs, MAXNDETECTIONS, run_detection_cascade,
								   pixels, nrows, ncols, ldim, scalefactor, stridefactor,
							 MAX(16, MIN_SIZE), MIN(128, MAX_SIZE));

		for(i=1; i<5; ++i)
		{
			cvResize(pyr[i-1], pyr[i], CV_INTER_LINEAR);

			pixels = (uint8_t*)pyr[i];
			nrows = height/(1<<i);	//1/(2^n)
			ncols = width/(1<<i);
			ldim = ncols;	//pyr[i]->widthStep;

			nd = find_objects(&rs[ndetections], &cs[ndetections], &ss[ndetections],
							  &qs[ndetections], MAXNDETECTIONS-ndetections,
					 run_detection_cascade, pixels, nrows, ncols, ldim, scalefactor,
					 stridefactor, MAX(64, MIN_SIZE>>i), MIN(128, MAX_SIZE>>i));

			for(j=ndetections; j<ndetections+nd; ++j)
			{
				rs[j] = (1<<i)*rs[j];
				cs[j] = (1<<i)*cs[j];
				ss[j] = (1<<i)*ss[j];
			}

			ndetections = ndetections + nd;
		}
		for(i=0; i<5; ++i){
			if(pyr[i]){
				free(pyr[i]);
				pyr[i]=NULL;
			}
		}
#endif
	}
	else
	{
		pixels = frame;//(uint8_t*)gray;
		nrows = height;
		ncols = width;
		ldim = width;	//gray->widthStep;
		printf("ldim=%d\n", ldim);
		//
		ndetections = find_objects(rs, cs, ss, qs, MAXNDETECTIONS, run_detection_cascade,
								   pixels, nrows, ncols, ldim, scalefactor,
							 stridefactor, MIN_SIZE, MIN(nrows, ncols));
	}

	ndetections = cluster_detections(rs, cs, ss, qs, ndetections);

#if 0
	t = getticks() - t;

	// if the flag is set, draw each detection

	if(draw)
		for(i=0; i<ndetections; ++i)
			if(qs[i]>=qthreshold) // check the confidence threshold
				cvCircle(frame, cvPoint(cs[i], rs[i]), ss[i]/2, CV_RGB(255, 0, 0), 4, 8, 0); // we draw circles here since height-to-width ratio of the detected face regions is 1.0f

	// if the flag is set, print the results to standard output
	if(print)
	{
		//
		for(i=0; i<ndetections; ++i)
			if(qs[i]>=qthreshold) // check the confidence threshold
				printf("%d %d %d %f\n", (int)rs[i], (int)cs[i], (int)ss[i], qs[i]);

		//
		printf("# %f\n", 1000.0f*t); // use '#' to ignore this line when parsing the output of the program
	}
#endif

	ndetections = (ndetections > maxdetect)?maxdetect:ndetections;
	for(i = 0, j=0; i < ndetections; i++)
		if(qs[i]>=qthreshold){
			frs[j]=rs[i];
			fcs[j]=cs[i];
			fss[j++]=ss[i];
	}
	return j;
}

