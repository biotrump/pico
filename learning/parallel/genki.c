#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <fcntl.h>
#include <errno.h>
#include <dirent.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "omp.h"
#include <cv.h>
#include <highgui.h>

#include "genki.h"

#define MAX_GENKI_RAND_BBOX		(7)
#define MAX_WIN_SIZE			(192.0)

int TotalObjListFiles;
int MaxObjFileNameLen;
void *pGenkiImgList;
void *pGenkiLabelList;
void *plistBuf;
/*
 * 0.0-1.0
 */
inline float rand_real(void)
{
	return rand()/(float)RAND_MAX;
}

/*
 * rid : raw intensity data
 *
 * path : the destination rid file and verification .bmp file for a face bounding box.
 * src : the grey level image of the genki face image.
 * roi : the newly cropped image ROI to the src image. The newly cropped box is a little bigger
 * 		around the face ROI.
 * ratio : the max size is MAX_WIN_SIZE pixel. If the column or height is greater than 192,
 * 		downsampling it to MAX_WIN_SIZE.
 */
int saveObjAsRID(char *path, IplImage* src,  CvRect roi, float ratio, GENKI_FACE_BBOX box)
{
	//assert(path && src);
	int fd = open(path, O_RDWR|O_CREAT|O_TRUNC, 0664 );
	_pr_debug("%s:roi(x,y,w,h)=(%d,%d,%d,%d)\n",__func__, roi.x, roi.y,
		   roi.width, roi.height);
	_pr_debug("src(w,h,s)=(%d,%d,%d)\n", src->width, src->height, src->imageSize);

	if(fd > 0){
		// Must have dimensions of output image
		IplImage* cropped = cvCreateImage( cvSize(roi.width,roi.height), src->depth, src->nChannels );
		IplImage *dst = NULL;
		if(cropped){
			// Say what the source region is
			cvSetImageROI( src, roi );
			// Do the copy
			cvCopy( src, cropped, NULL );
			if(ratio<1.0){//resize the cropped image
				_pr_debug("src:%d x %d : %d bytes\n", src->width, src->height,
						src->imageSize);
				_pr_debug("crop:%d x %d : %d bytes\n", cropped->width, cropped->height,
						cropped->imageSize);
				CvSize dst_cvsize = cvSize(cropped->width * ratio, cropped->height * ratio);
				dst = cvCreateImage( dst_cvsize, cropped->depth, cropped->nChannels);
				cvResize(cropped, dst, CV_INTER_LANCZOS4);
				_pr_debug("dst:%d x %d : %d bytes\n", dst->width, dst->height,
					dst->imageSize);
			}
			cvResetImageROI( src );
			CvPoint VertexOne=cvPoint(box.center_col-box.diameter/2,
									  box.center_row-box.diameter/2);
			CvPoint VertexThree=cvPoint(box.center_col+box.diameter/2,
									  box.center_row+box.diameter/2);
			if(dst){
				strcpy(path + strlen(path)-4, "_r.bmp");
				cvRectangle(dst, VertexOne, VertexThree, CV_RGB(255,0, 0), 2, CV_AA, 0);
				cvSaveImage(path , dst, 0);//rid's bmp file for reference.
				cvReleaseImage(&dst);
			}else{
				strcpy(path + strlen(path)-3, "bmp");
				cvRectangle(cropped, VertexOne, VertexThree, CV_RGB(255,0,0), 2, CV_AA, 0);
				cvSaveImage(path , cropped, 0);//rid's bmp file for reference.
			}
			write(fd, &roi.width,sizeof(int));	//4 bytes w
			write(fd, &roi.height,sizeof(int));	//4 bytes h
			write(fd, cropped->imageData, cropped->imageSize);	//body
			_pr_debug("%d x %d : %d bytes\n", roi.width, roi.height, cropped->imageSize);
			cvReleaseImage(&cropped);

		}else{
			printf("[%d,%d,%d,%d] crop fails to create.\n", roi.width,roi.height,
				src->depth, src->nChannels);
		}
		close(fd);
	}else{
		printf("[%s] fails to open:%d\n", path,errno);
	}
	return errno;
}

/*
 * write 7 random boxes to list.txt for each rid.
 *  1 face0.rid
 *  2     7	==> 7 random boxs
 *  3     54.128180 54.193473 61.633318
 *  4     50.617332 47.980850 62.638414
 *  5     54.102598 52.237606 61.709570
 *  6     53.919454 51.292225 71.041327
 *  7     52.621050 51.479368 70.888629
 *  8     50.815017 53.567457 62.003842
 *  9     51.473735 49.503986 67.789544
 *
 */
void export(char *dstfolder, IplImage *image, GENKI_FACE_BBOX bbox, int id, FILE *fList)
{
	assert(image);
	CvRect roi;
	int nrows = image->height;
	int ncols = image->width;
	int r, c, s;
	r = bbox.center_row;
	c = bbox.center_col;
	s = bbox.diameter;

	/* The new cropped area is bounded by +-0.75diameter around box center (r,c).
	 * The new cropped area is bigger than the original face bounding box, +-0.5diameter,
	 * around box center (r,c)
	 */
	int r0 = max((r - 0.75*s), 0);
	int r1 = min(r + 0.75*s, nrows);
	int c0 = max((c - 0.75*s), 0);
	int c1 = min(c + 0.75*s, ncols);

	//the new cropped image containing the original face bounding box.
	//the new dim of the cropped image
	nrows = r1- r0 + 1;
	ncols = c1- c0 + 1;

	//the new bounding box center in the new cropped image
	r = r - r0;
	c = c - c0;

	//resize, if needed
	int wsize = max(nrows, ncols);
	float ratio = MAX_WIN_SIZE/wsize;

	if (ratio<1.0){//resize the pic because it's > MAX_WIN_SIZE
		//im = numpy.asarray( Image.fromarray(im).resize((int(ratio*ncols), int(ratio*nrows))) )
		_pr_debug("++++ratio=%f\n",ratio);
		r = ratio*r;
		c = ratio*c;
		s = ratio*s;
	}
#if 1
	char (*listBuf)[(MAX_GENKI_RAND_BBOX+1)*MAX_OBJ_FILENAME_LEN]=plistBuf;
	if(listBuf){
		char *strBuf= (char *)(listBuf+id);
		_pr_debug("id=%d\n", id);
		memset(strBuf, 0, (MAX_GENKI_RAND_BBOX+1)*MAX_OBJ_FILENAME_LEN);
		sprintf(strBuf, "face%d.rid\n",id);
		strBuf += strlen(strBuf);
		//creating 7 randomized face bounding boxes from the original bounding box
		int nrands = MAX_GENKI_RAND_BBOX;
		//write '7' to list.txt
		sprintf(strBuf, "\t%d\tr=%.3f\tr=%d\tc=%d\ts=%d\n",nrands,
				ratio,r,c, s);
		strBuf += strlen(strBuf);
		int i;
		for(i= 0; i < nrands; i ++){
			float stmp, rtmp, ctmp;
			//uniformly randomize size (diameter of the bounding box) ratio 0.9-1.1, 10% random
			//stmp = s*random.uniform(0.9, 1.1)
			//stmp = s * (0.9 + 0.2 * cvRandReal(&rng));
			stmp = s * (0.9 + 0.2 * rand_real());
			//uniformly randomize row and column, ratio +-5% random
			//rtmp = r + s*random.uniform(-0.05, 0.05)
			//rtmp = r + s*(0.05 - 0.1 *cvRandReal(&rng));
			rtmp = r + s*(0.05 - 0.1 *rand_real());
			//ctmp = c + s*random.uniform(-0.05, 0.05)
			//ctmp = c + s*(0.05 - 0.1 *cvRandReal(&rng));
			ctmp = c + s*(0.05 - 0.1 *rand_real());
			//write the randomized row, column and size (diameter of bounding box)
			sprintf(strBuf, "\t%f %f %f\n",rtmp, ctmp, stmp);
			strBuf += strlen(strBuf);
		}
	}
#else
	#pragma omp critical
	{
		fprintf(fList, "face%d.rid\n",id);

		//creating 7 randomized face bounding boxes from the original bounding box
		int nrands = MAX_GENKI_RAND_BBOX;
		//write '7' to list.txt
		fprintf(fList, "\t%d\tr=%.3f\n",nrands, ratio);
		/*fprintf(fList, "\t[r0=%d\tr1=%d\tc0=%d\tc1=%d]\t[nrows=%d\tncols=%d]\n",
				r0, 	r1, c0, c1, nrows, ncols);
		fprintf(fList, "\t[%d\t%d\t%d]=>\t[r=%d\tc=%d\ts=%d]\n",
				bbox.center_row, bbox.center_col, bbox.diameter, r, c, s);
		fprintf(fList, "\twsize=%d\tratio=%f\n",wsize, ratio);*/
		int i;
		for(i= 0; i < nrands; i ++){
			float stmp, rtmp, ctmp;
			//uniformly randomize size (diameter of the bounding box) ratio 0.9-1.1, 10% random
			//stmp = s*random.uniform(0.9, 1.1)
			//stmp = s * (0.9 + 0.2 * cvRandReal(&rng));
			stmp = s * (0.9 + 0.2 * rand_real());
			//uniformly randomize row and column, ratio +-5% random
			//rtmp = r + s*random.uniform(-0.05, 0.05)
			//rtmp = r + s*(0.05 - 0.1 *cvRandReal(&rng));
			rtmp = r + s*(0.05 - 0.1 *rand_real());
			//ctmp = c + s*random.uniform(-0.05, 0.05)
			//ctmp = c + s*(0.05 - 0.1 *cvRandReal(&rng));
			ctmp = c + s*(0.05 - 0.1 *rand_real());
			//write the randomized row, column and size (diameter of bounding box)
			fprintf(fList, "\t%f %f %f\n",rtmp, ctmp, stmp);
		}
		fputc('\n',fList);
	}
#endif
	char ridFileName[MAX_FILE_PATH_SIZE];
	roi = cvRect(c0, r0, ncols-1, nrows-1);
	snprintf(ridFileName, MAX_FILE_PATH_SIZE, "%s/face%d.rid",dstfolder, id);
	//printf("ridFileName=%s\n", ridFileName);
	bbox.center_row=r;
	bbox.center_col=c;
	bbox.diameter=s;
	saveObjAsRID(ridFileName, image, roi, ratio, bbox);
}

void exportmirrored(char *dstfolder, IplImage *image, GENKI_FACE_BBOX bbox, int id, FILE *fList)
{
	/*
	 * exploit mirror symmetry of the face
	 *
	 * flip image, horizontal , so the dim is not changed!
	 */
	//im = numpy.asarray(ImageOps.mirror(Image.fromarray(im)))
	cvFlip(image, image, 1);
	// flip column, so the face bounding box is mirrored.
	// the center of the bounding box is shifted in x-axis after the mirror, too.
	//c = im.shape[1] - c
	bbox.center_col = image->width - bbox.center_col ;
	// export
	export(dstfolder, image, bbox, id, fList);
}

/*
 * read the object file list and store the list
 */
int readGenkiImgFileList(char *path)
{
	FILE *f;
	char t[1024], filename[MAX_FILE_PATH_SIZE];
	int ret=-1, len=1;
	assert(path);
	snprintf(filename, MAX_FILE_PATH_SIZE, "%s%s", path, CV_DATASET_GENKI_IMAGES);

	//counting the lines in the file.
	if((f=fopen(filename, "rt"))!=NULL){
		while(fgets(t,1024,f) != NULL){
		TotalObjListFiles++;
		if(strnlen(t, 1024) > len) len = strnlen(t, 1024);
	}
	MaxObjFileNameLen=len+1;//plus null
	rewind(f);
	printf("total %d lines, max strlen=%d\n", TotalObjListFiles, MaxObjFileNameLen);

	if(TotalObjListFiles && MaxObjFileNameLen){
	  int i = 0;
	  pGenkiImgList = calloc(TotalObjListFiles, MaxObjFileNameLen);
	  char (*p)[MaxObjFileNameLen] = pGenkiImgList;
	  while((i < TotalObjListFiles) &&
			fgets(p[i],MAX_OBJ_FILENAME_LEN,f) != NULL){
		//remove \n or \r
		int l=strnlen(p[i], 1024);
		if(l){
		  if( (p[i][l-1] == '\n') || (p[i][l-1] == '\r'))
			p[i][l-1]=0;
		  _pr_debug("%d:%s\n",i, p[i]);
		}
	  	i++;
	  }
	  printf("total %d files in the list:%s\n",i, p[i-1]);
	}
	fclose(f);
	ret = 0;
  }
  return ret;
}

/*
 * GENKI/GENKI-R2009a/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt
 * Each line is a bounding box of a face in the same line number of GENKI-SZSL_Images.txt.
 * For example: the first line of GENKI-SZSL_Images.txt is "file0000000000003987.jpg",
 * so the first line of GENKI-SZSL_labels.txt is the bounding box of the file file0000000000003987.jpg.
 *
 * read the face bounding box (x,y,diameter) mapped to each file listed in GENKI-SZSL_Images.txt
 */
int readGenkiLabelList(char *path)
{
	FILE *f;
	char str[MAX_FILE_PATH_SIZE], filename[MAX_FILE_PATH_SIZE];
	int ret=-1, bc=0;
	assert(path);
	snprintf(filename, MAX_FILE_PATH_SIZE, "%s%s", path, CV_DATASET_GENKI_LABELS);

	if((f=fopen(filename, "rt"))!=NULL){
		while(fgets(str,MAX_FILE_PATH_SIZE,f) != NULL){
			bc++;
		}
		rewind(f);
		printf("total %d box lines\n", bc);
		if(bc){
			int i = 0;
			pGenkiLabelList = calloc(bc, sizeof(GENKI_FACE_BBOX));
			PGENKI_FACE_BBOX p = (PGENKI_FACE_BBOX)pGenkiLabelList;
			while(	fscanf(f, "%d %d %d", &(p[i].center_col), &(p[i].center_row),
						&(p[i].diameter)) != EOF){
				_pr_debug("%d:(x,y,d)=(%d, %d, %d)\n",i, p[i].center_col,p[i].center_row,p[i].diameter );
				i++;
			}
			printf("total %d files in the list:(%d,%d,%d)\n",bc,
				p[i-1].center_col, p[i-1].center_row, p[i-1].diameter);
		}
		fclose(f);
		ret = 0;
	}else{
		printf("%s opening failure\n", filename);
	}
	return ret;
}

void memdump(void)
{
	struct sysinfo info;
	sysinfo( &info );
	printf("mem_unit=%lu\n", (size_t)info.mem_unit);
	printf("total RAM=%lu\n", (size_t)info.totalram * (size_t)info.mem_unit);
	printf("total Available RAM=%lu\n", (size_t)info.freeram * (size_t)info.mem_unit);
}

int genkiFace(char *facepath,char *riddst)
{
	unsigned long ulUsed=0;
	int curObjIndex=0;
	int totalNonObjFiles=0;

	memdump();
	if( facepath && riddst){
		char listFileName[MAX_FILE_PATH_SIZE];
		struct stat st = {0};
		if (stat(riddst, &st) == -1) {
			mkdir(riddst, 0766);
			printf("mkdir :%s\n", riddst);
		}
		//char *srcfolder,*dstfolder;
		readGenkiImgFileList(facepath);
		readGenkiLabelList(facepath);

		snprintf(listFileName, MAX_FILE_PATH_SIZE, "%s/list.txt",riddst);
		FILE *fList = fopen(listFileName, "w" );
		//total files including original and mirrored obj files.
		plistBuf = calloc(TotalObjListFiles*2, (MAX_GENKI_RAND_BBOX+1)*MAX_OBJ_FILENAME_LEN);
		if(plistBuf){
			if(fList > 0){
				fprintf(fList, "%d\n",TotalObjListFiles*2);
			#pragma omp parallel
			{
				int id = omp_get_thread_num();//local to this thread
				int counts=0,i=0;
				int percent=0;
				int cont=1;	//local to this thread stack
				char (*list)[MaxObjFileNameLen] = pGenkiImgList;
				char imgFileName[MAX_FILE_PATH_SIZE];
				PGENKI_FACE_BBOX fb=NULL;

				while(cont){
					memset(imgFileName,0,MAX_FILE_PATH_SIZE);
					#pragma omp critical
					{//get an id to the file list, so the file can be opened later.
						if(curObjIndex < TotalObjListFiles){
							snprintf(imgFileName, MAX_FILE_PATH_SIZE, "%s/files/%s",facepath, list[curObjIndex]);
							i=curObjIndex++;
							percent = curObjIndex*1000.0 / TotalObjListFiles;
							_pr_debug(">>>CR:%d:%d:%s\n", id, i, imgFileName);
						}
						if(curObjIndex >= TotalObjListFiles){
							printf("%d is done\n", id);
							cont=0;
						}
					}
					if(imgFileName[0]){//load the image file
						IplImage *image=NULL;
						fb = (PGENKI_FACE_BBOX)pGenkiLabelList;
						if(0 == (percent % 10)) printf("%3d%%\n", percent/10);
						_pr_debug("i->%d\n", i);
						_pr_debug("<<<TID[%d]:%d:(%d,%d,%d):%s\n", id, i,
							fb[i].center_col, fb[i].center_row, fb[i].diameter, imgFileName);
						image= cvLoadImage(imgFileName, CV_LOAD_IMAGE_GRAYSCALE);
						if(image){
							_pr_debug(">>>i->%d\n", i);
							_pr_debug("depth:%d, channel=%d, (%d,%d)\n", image->depth, image->nChannels,
															image->width, image->height);
							export(riddst, image, fb[i] , i*2, fList);
							exportmirrored(riddst, image, fb[i] , i*2+1, fList);
							cvReleaseImage(&image);
						}else{
							printf("%s not found\n", imgFileName);
							cont=0;
						}
					}else{
						printf("TID:%d, null image\n", id);
					}
				}
			}
			int j;
			char (*listBuf)[(MAX_GENKI_RAND_BBOX+1)*MAX_OBJ_FILENAME_LEN]=plistBuf;
			/* TotalObjListFiles * 2 : original and mirrored obj*/
			for(j = 0; j < TotalObjListFiles*2; j ++){
				fputs(listBuf[j], fList);
				_pr_debug("%s",listBuf[j]);
			}
			fclose(fList);

			}else{
				printf("[%s]opening failure : %d\n", listFileName, errno);
			}
		if(plistBuf)
			free(plistBuf);
		}else{
			printf("listBuf allocation error\n");
		}
		if(pGenkiImgList){
			free(pGenkiImgList);
			pGenkiImgList=NULL;
		}
		if(pGenkiLabelList){
			free(pGenkiLabelList);
			pGenkiLabelList=NULL;
		}
	}else{
		printf("facepath or riddst is empty...\n");
	}

	memdump();
}
