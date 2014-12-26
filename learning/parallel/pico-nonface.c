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

int TotalNonObjListFiles;
int MaxNonObjFileNameLen;
void *pNonObjFileList;

/*
 * rid : raw intensity data
 *
 * path : the destination rid file and verification .bmp file for a face bounding box.
 * src : the grey level image of the genki face image.
 * roi : the newly cropped image ROI to the src image. The newly cropped box is a little bigger
 * 		around the face ROI.
 * ratio : the max size is MAX_WIN_SIZE pixel. If the column or height is greater than 192,
 * 		downsampling it to MAX_WIN_SIZE.
 * Some jpg has EXIF meta data which has orientation value, but it's ignored by cvLoadImage.
 */
int saveNonObjAsRID(char *path, IplImage* src, int i)
{
	assert(path && src);
	//printf("%s:%s\n",__func__, path);
	//printf("src(w,h,s)=(%d,%d,%d)\n", src->width, src->height, src->imageSize);
	char (*list)[MaxNonObjFileNameLen] = pNonObjFileList;
	char nonObjRIDName[MAX_FILE_PATH_SIZE];
	// Must have dimensions of output image
	snprintf(nonObjRIDName, MAX_FILE_PATH_SIZE, "%s/%05d-%s.rid",
				 path, i, list[i]);
	//printf("nonObjRIDName=%s\n", nonObjRIDName);
	int fd = open(nonObjRIDName, O_RDWR|O_CREAT|O_TRUNC, 0664 );
	if(fd > 0){
		strcpy(nonObjRIDName + strlen(nonObjRIDName)-7, "bmp");
		//printf("+++nonObjRIDName=%s\n", nonObjRIDName);
		cvSaveImage(nonObjRIDName , src, 0);	//rid's bmp file for reference.
		write(fd, &src->width, sizeof(int));	//4 bytes w
		write(fd, &src->height, sizeof(int));	//4 bytes h
		write(fd, src->imageData, src->imageSize);	//body
		//printf("%d x %d : %d bytes\n", roi.width, roi.height, cropped->imageSize);
		close(fd);
	}else{
		printf("[%s] fails to open:%d\n", path, errno);
	}
	return errno;
}


/*
 * read the non object file and create a name list
 */
int scanNonObjFile(char *path, FILE *fList)
{
	assert(path);
	int ret=0;
	int i=0, n=0;
	int maxl=0;
	DIR *d;
	struct dirent *dir;
	if (d = opendir(path))
	{
		while ((dir = readdir(d)) != NULL)
		{
			if(strncmp(".", dir->d_name, strlen(".")) &&
				strncmp("..", dir->d_name, strlen(".."))){
				n++;
				if(strlen(dir->d_name) > maxl) maxl=strlen(dir->d_name);
			}
			_pr_debug("%d: %s\n", n,  dir->d_name);
		}
		TotalNonObjListFiles=n;
		rewinddir(d);
		printf("nonobj file name max leng:%d\n",maxl);
		maxl++;
		MaxNonObjFileNameLen = maxl;
		pNonObjFileList=calloc(n, maxl);
		char (*nonObj)[maxl] = pNonObjFileList;
		char nonObjRIDName[20+maxl];
		if(pNonObjFileList){
			fprintf(fList, "%d\n",TotalNonObjListFiles);
			while ( (i < n) && ((dir = readdir(d)) != NULL))
			{
				if(strncmp(".", dir->d_name, strlen(".")) &&
					strncmp("..", dir->d_name, strlen(".."))){
					strncpy(nonObj[i], dir->d_name, strlen(dir->d_name));
					_pr_debug("%d,len=%zd: %s\n", i, strlen(nonObj[i]), nonObj[i]);
					//write entry 00000-IMG_2121.JPG.rid to list.txt.
					memset(nonObjRIDName, 0, 20+maxl);
					snprintf(nonObjRIDName, 20+maxl, "%05d-%s.rid\n", i, nonObj[i]);
					fputs(nonObjRIDName, fList);
					_pr_debug(">>:%s\n", nonObjRIDName);
					i ++;
				}
			}
		}
		closedir(d);
	}else{
		printf("%s opening failure\n", path);
		ret = -1;
	}
  return ret;
}

int picoNonFace(char *nonfacepath, char *nonObjrid)
{
	int i, curNonObjIndex=0;

	memdump();
	if( nonfacepath && nonObjrid){
		char listFileName[MAX_FILE_PATH_SIZE];
		struct stat st = {0};
		if (stat(nonObjrid, &st) == -1) {
			mkdir(nonObjrid, 0766);
			printf("mkdir :%s\n", nonObjrid);
		}
		snprintf(listFileName, MAX_FILE_PATH_SIZE, "%s/list.txt", nonObjrid);
		FILE *fList = fopen(listFileName, "w+" );
		if(fList > 0){
			scanNonObjFile(nonfacepath, fList);
			printf("starting processing nonObj image by omp...\n");
			#pragma omp parallel
			{
				int id = omp_get_thread_num();//local to this thread
				int i=0;
				int percent=0;
				int cont=1;	//local to this thread stack
				char (*list)[MaxNonObjFileNameLen] = pNonObjFileList;
				char nonObjFilePath[MAX_FILE_PATH_SIZE];

				while(cont){
					memset(nonObjFilePath,0,MAX_FILE_PATH_SIZE);
					#pragma omp critical
					{//get an id to the file list, so the file can be opened later.
						if(curNonObjIndex < TotalNonObjListFiles){
							snprintf(nonObjFilePath, MAX_FILE_PATH_SIZE, "%s/%s",
									nonfacepath, list[curNonObjIndex]);
							_pr_debug(">>>CR:tid[%d]:%d:%s\n", id, curNonObjIndex, nonObjFilePath);
							i = curNonObjIndex++;
							percent = curNonObjIndex*1000.0 / TotalNonObjListFiles;
						}
						if(curNonObjIndex >= TotalNonObjListFiles){
							printf("%d is done\n", id);
							cont=0;
						}
					}
					if(nonObjFilePath[0]){//load the image file
						_pr_debug("i->%d\n", i);
						_pr_debug("<<<TID[%d]:%s\n", id, nonObjFilePath);
						if(0 == (percent % 10)) printf("%3d%%\n", percent/10);
						IplImage *image=cvLoadImage(nonObjFilePath, CV_LOAD_IMAGE_GRAYSCALE);
						if(image){
							saveNonObjAsRID(nonObjrid, image, i);
							cvReleaseImage(&image);
						}else{
							printf("%s cvLoadImage failure\n", nonObjFilePath);
							cont=0;
						}
					}
				}
			}
			fclose(fList);
		}else{
			printf("[%s]opening failure : %d\n", listFileName, errno);
		}
		if(pNonObjFileList){
			free(pNonObjFileList);
			pNonObjFileList=NULL;
		}
	}else{
		printf("nonfacepath or nonObjrid is empty!!!\n");
	}

	memdump();
}
