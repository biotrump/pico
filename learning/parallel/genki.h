#ifndef	_H_GENKI_HEADER_H
#define	_H_GENKI_HEADER_H

#define	MAX_OBJ_FILENAME_LEN	(81)
#define	MAX_OBJ_FILES		(1000000)
#define	MAX_FILE_PATH_SIZE	(1024)

#define	CV_DATASET_GENKI_IMAGES		"/Subsets/GENKI-SZSL/GENKI-SZSL_Images.txt"
#define	CV_DATASET_GENKI_LABELS	   	"/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt"

#define	PICO_INIT_TSR		(1000)
#define	PICO_INIT_TSC		(1001)
#define	PICO_INIT_STAGAES	(1002)

#define max( a, b ) ((a) > (b) ? (a) : (b) )
#define min( a, b ) ((a) < (b) ? (a) : (b) )
#define MAX_WIN_SIZE	(192.0)

#define _pr_debug(fmt, arg...) \
	if (Verbose) fprintf(stderr, fmt, ##arg)

typedef struct _genki_bounding_box{
int	center_col;	//X of the box center
int center_row;	//Y of the box center
int diameter;	//diameter of the box
} GENKI_FACE_BBOX, *PGENKI_FACE_BBOX;

typedef struct _lrn_param{
	char *src_path;
	char *faces_path;
	char *nonfaces_path;
	int maxnstages;
	float targetfpr;
	int tdepths;
	float minstagetpr;
	float maxstagefpr;
	int maxnumtreesperstage;
	char *dst_path;
} LEARN_PARAM, *PLEARN_PARAM;

typedef struct _lrn_param_tbl{
	PLEARN_PARAM plrn;
	int totalParams;
} LEARN_PARAM_TBL, *PLEARN_PARAM_TBL;

extern int genkiFace(char *facepath, char *riddst);
extern int picoNonFace(char *nonfacepath, char *nonObjrid);
extern void memdump(void);

extern int TotalListFiles;
extern int MaxFileNameLen;
extern void *pGenkiImgList;
extern void *pGenkiLabelList;
extern int	Verbose;
#endif
