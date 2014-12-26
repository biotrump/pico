//#include <cv.h>
//#include <highgui.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "omp.h"

#if 0
           struct sysinfo {
               long uptime;             /* Seconds since boot */
               unsigned long loads[3];  /* 1, 5, and 15 minute load averages */
               unsigned long totalram;  /* Total usable main memory size */
               unsigned long freeram;   /* Available memory size */
               unsigned long sharedram; /* Amount of shared memory */
               unsigned long bufferram; /* Memory used by buffers */
               unsigned long totalswap; /* Total swap space size */
               unsigned long freeswap;  /* swap space still available */
               unsigned short procs;    /* Number of current processes */
               unsigned long totalhigh; /* Total high memory size */
               unsigned long freehigh;  /* Available high memory size */
               unsigned int mem_unit;   /* Memory unit size in bytes */
               char _f[20-2*sizeof(long)-sizeof(int)]; /* Padding for libc5 */
           };
#endif

#define	MAX_OBJ_FILENAME_LEN	(81)
#define	MAX_OBJ_FILES		(1000000)
#define	MAX_FILE_PATH_SIZE	(4096)

#define	CV_DATASET_GENKI_IMAGES		"/Subsets/GENKI-SZSL/GENKI-SZSL_Images.txt"
#define	CV_DATASET_GENKI_LABELS	   	"/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt"

typedef struct _genki_bounding_box{
int	center_col;	//X of the box center
int center_row;	//Y of the box center
int diameter;	//diameter of the box
} GENKI_FACE_BBOX, *PGENKI_FACE_BBOX;

int TotalListFiles;
int MaxFileNameLen;
void *pGenkiImgList;
void *pGenkiLabelList;

/*
 * raw intensity data
 *
 * im is an array/list from a cropped image around the face bounding box:
 * im dim : row * col
 * rid file : 4 bytes width, 4 bytes height , image data (row-wise)
 */
int saveasrid(char *path, int w, int h)
{
  int fd = open(path, O_RDWR|O_CREAT,0664 );
  if(fd > 0){
	write(fd, &w, sizeof(w));
	write(fd, &h, sizeof(h));

	close(fd);
  }
  return errno;
}

/*
 * read the object file list and store the list
 */
int fnGenkiImgFileList(char *path)
{
  FILE *f;
  char t[1024], filename[MAX_FILE_PATH_SIZE];
  int ret=-1, len=1;
  assert(path);
  snprintf(filename, MAX_FILE_PATH_SIZE, "%s%s", path, CV_DATASET_GENKI_IMAGES);

  if((f=fopen(filename, "rt"))!=NULL){
	while(fgets(t,1024,f) != NULL){
	  TotalListFiles++;
	  if(strnlen(t, 1024) > len) len = strnlen(t, 1024);
	}
	MaxFileNameLen=len+1;//plus null
	rewind(f);
	printf("total %d lines, max strlen=%d\n", TotalListFiles, MaxFileNameLen);
	if(TotalListFiles && MaxFileNameLen){
	  int i = 0;
	  pGenkiImgList = calloc(TotalListFiles, MaxFileNameLen);
	  char (*p)[MaxFileNameLen] = pGenkiImgList;
	  while((i < TotalListFiles) &&
			fgets(p[i],MAX_OBJ_FILENAME_LEN,f) != NULL){
		//remove \n or \r
		int l=strnlen(p[i], 1024);
		if(l){
		  if( (p[i][l-1] == '\n') || (p[i][l-1] == '\r'))
			p[i][l-1]=0;
		  //printf("%d:%s\n",i, p[i]);
		}
	  	i++;
	  }
	  printf("%s\n",p[0]);
	  printf("total %d files in the list\n",i);
	}
	fclose(f);
	ret = 0;
  }
  char (*p)[MaxFileNameLen] = pGenkiImgList;
  printf("%p %s\n",pGenkiImgList, p[0]);
  return ret;
}

/*
 * read the bouding box (x,y,diameter) mapped to each file listed in images.txt
 */
int fnGenkiList(char *path)
{
  FILE *f;
  char t[1024], filename[MAX_FILE_PATH_SIZE];
  int ret=-1, bc=0;
  assert(path);
  snprintf(filename, MAX_FILE_PATH_SIZE, "%s%s", path, CV_DATASET_GENKI_LABELS);

  if((f=fopen(filename, "rt"))!=NULL){
	while(fgets(t,1024,f) != NULL){
	  bc++;
	}
	rewind(f);
	printf("total %d box lines\n", bc);
	if(bc){
	  int i = 0;
	  pGenkiLabelList = calloc(bc, sizeof(GENKI_FACE_BBOX));
	  PGENKI_FACE_BBOX p = (PGENKI_FACE_BBOX)pGenkiImgList;
	  while(	fscanf(f, "%d %d %d", &(p[i].center_col), &(p[i].center_row),
				   &(p[i].diameter)) != EOF){
		  //printf("%d:(x,y,d)=(%d, %d, %d)\n",i, p[i].center_col,p[i].center_row,p[i].diameter );
		  i++;
		}
	}
	printf("total %d files in the list\n",bc);
	fclose(f);
	ret = 0;
  }
  return ret;
}

int meminfo()
{
  struct sysinfo info;
  sysinfo( &info );
  printf("mem_unit=%lu\n", (size_t)info.mem_unit);
  printf("total RAM=%lu\n", (size_t)info.totalram * (size_t)info.mem_unit);
  printf("total Available RAM=%lu\n", (size_t)info.freeram * (size_t)info.mem_unit);
}

int main(int argc, char **argv)
{
  unsigned long ulUsed=0;
  static int curObjIndex=0;
  int totalNonObjFiles=0;
  /*CvRNG rng;
  rng = cvRNG(cvGetTickCount());*/
  //read the file list lines to curObjIndex
  meminfo();

  if( (argc >= 2) && argv[1]){
	fnGenkiImgFileList(argv[1]);
	fnGenkiList(argv[1]);

  char (*l)[MaxFileNameLen] = pGenkiImgList;
  //PGENKI_FACE_BBOX ll = (PGENKI_FACE_BBOX)pGenkiImgList;
  printf("%p %p %d, total entries %d, %s\n",pGenkiImgList, l, sizeof(pGenkiImgList) ,TotalListFiles, l[0]);
  //printf("(%d,%d,%d)\n",ll[0].center_col, ll[0].center_row, ll[0].diameter);


  }

  exit(1);
  meminfo();
  //allocating memory till 1G left!
  if(pGenkiImgList && pGenkiLabelList){
	printf(">>>\n");
	#pragma omp parallel
	{
	  int id = omp_get_thread_num();//local to this thread
	  int counts=0,i=0;
	  int cont=1;	//local to this thread stack
	  char (*list)[MaxFileNameLen] = pGenkiImgList;
	  char myFileName[MAX_FILE_PATH_SIZE];
	  memset(myFileName,0,MAX_FILE_PATH_SIZE);

	  while(cont){
		#pragma omp critical
		{//get an id to the file list, so the file can be opened later.
		  if(curObjIndex < TotalListFiles){
			  //printf("%d : curObjIndex=%d, [%s]\n",id, curObjIndex, list[curObjIndex]);
			  //strncpy(myFileName, list[curObjIndex], MaxFileNameLen);
			  snprintf(myFileName, MAX_FILE_PATH_SIZE, "%s/files/%s",argv[1], list[curObjIndex]);
			  i=curObjIndex++;
		  }else{
			  printf("%d is done\n", id);
			  cont=0;
		  }
		}
		if(cont){
		  int fd;
		  printf("[%d][%d]:%s is processing...\n",id, i, myFileName);
		  if( (fd = open(myFileName, O_RDONLY)) > 0){
			printf("%s is opening\n", myFileName);
			counts++;
			//sleep(random()%1);
			close(fd);
		  }else{
			printf("%s opening failure\n", myFileName, errno);
			exit(1);
		  }
		}
	  }
	  printf("<<<%d exits, total files [%d] processed.\n", id, counts);
	}
	printf("total %d exit...\n", curObjIndex);
  }

  if(pGenkiImgList){
	free(pGenkiImgList);
	pGenkiImgList=NULL;
  }
  if(pGenkiLabelList){
	free(pGenkiLabelList);
	pGenkiLabelList=NULL;
  }
  meminfo();

}