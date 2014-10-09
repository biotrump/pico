/*
 *	Copyright (c) 2013, Nenad Markus
 *	All rights reserved.
 *
 *	This is an implementation of the algorithm described in the following paper:
 *		N. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer,
 *		A method for object detection based on pixel intensity comparisons,
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

#include <stdio.h>
#include <time.h> //!!!don't use "sys/time.h" for C
#include <sys/stat.h>
//#include <cv.h>
//#include <highgui.h>
#include <errno.h>
#include <fcntl.h>

#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "capture-v4l2.h"
//#define _ROTATION_INVARIANT_DETECTION_

/*

*/

#ifndef _ROTATION_INVARIANT_DETECTION_
#define _INLINE_BINTEST_
#endif

#include "../../picort.c"

/*
	object detection parameters
*/
#define	 FRAME_WIDTH	(640)
#define FRAME_HEIGHT	(480)
#define FRAME_FPS		(25)

#ifndef QCUTOFF
#define QCUTOFF 3.0f
#endif

#ifndef MINSIZE
#define MINSIZE 100
#endif

#ifndef SCALEFACTOR
#define SCALEFACTOR 1.2f
#endif

#ifndef STRIDEFACTOR
#define STRIDEFACTOR 0.1f
#endif

#ifndef MAXNDETECTIONS
#define MAXNDETECTIONS 10	//2048 is too many
#endif

#define	ROI_BORDERW	(4)
#define	ROI_WIDTH	(160)
#define	ROI_HEIGHT	(160)

/*

*/

int minsize = MINSIZE;

/* openCV cvCaptureFromCAM will lock the /dev/videoxxx, so the VIDIOC_S_FMT
can't work.
Set the capture format before openCV open the camera device.
*/
int preset_videofmt(int idx)
{
	struct v4l2_format fmt;
	struct v4l2_frmivalenum frmival;
	uint32_t fps;
	int camfd;
	char dev[20]={0};
	sprintf(dev,"/dev/video%d",idx);
	camfd = open(dev, O_RDWR);
	if (camfd == -1){
        perror("Opening video device");
        return -1;
	}
	print_caps(camfd);
	EnumVideoFMT(camfd);
	GetVideoFMT(camfd, &fmt);
//    fmt.fmt.pix.width = 640;
//    fmt.fmt.pix.height = 480;
	fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV; //V4L2_PIX_FMT_MJPEG;
//	fmt.fmt.pix.bytesperline = 1280;
//	fmt.fmt.pix.sizeimage = 614400;
//	fmt.fmt.pix.colorspace=0;
//	fmt.fmt.pix.field = V4L2_FIELD_NONE;
	SetVideoFMT(camfd, fmt);
	GetVideoFMT(camfd, &fmt);
	EnumFrameRate(camfd, V4L2_PIX_FMT_YUYV);

	memset(&frmival,0,sizeof(frmival));
    frmival.pixel_format = V4L2_PIX_FMT_YUYV;
    frmival.width = FRAME_WIDTH;
    frmival.height = FRAME_HEIGHT;
	fps = GetFPSParam(camfd, 9.0, &frmival);
	SetFPSParam(camfd, fps);

	close(camfd);
	//15fps
}

int CamParameter(int camfd)
{
	struct v4l2_format fmt;
	//disabling AWB.
	//printf("fd=%d\n",camfd);
//	print_caps(camfd);
	printf("AWB=%d\n",GetAutoWhiteBalance(camfd));
//	SetAutoWhiteBalance(camfd,0);
	printf("AWB=%d\n",GetAutoWhiteBalance(camfd));
//	EnumVideoFMT(camfd);

	//15fps
}

int postset_cam(int camfd)
{

	CamParameter(camfd);
	printf("fd=%d\n",camfd);
	GetAutoExposure(camfd);
	SetAutoExposure(camfd, /*V4L2_EXPOSURE_MANUAL ,*/ V4L2_EXPOSURE_APERTURE_PRIORITY  );
	SetAutoExposureAutoPriority(camfd,0);
	printf("AutoPriority=%d\n",GetAutoExposureAutoPriority(camfd));
//	PrintFrameInterval(camfd, frmival.pixel_format, frmival.width, frmival.height);
	SetManualExposure(camfd, 140);
	GetManualExposure(camfd);
}

IplImage* cvGetSubImage(IplImage *img1, CvRect roi)
{
	/* load image */
	//IplImage *img1 = cvLoadImage("elvita.jpg", 1);

	/* sets the Region of Interest
	   Note that the rectangle area has to be __INSIDE__ the image */
	//cvSetImageROI(img1, cvRect(10, 15, 150, 250));
	cvSetImageROI(img1,roi);
	/* create destination image
	   Note that cvGetSize will return the width and the height of ROI */
	IplImage *img2 = cvCreateImage(cvGetSize(img1),
		                           img1->depth,
		                           img1->nChannels);

	/* copy subimage */
	cvCopy(img1, img2, NULL);

	/* always reset the Region of Interest */
	cvResetImageROI(img1);

}

inline void cvFreeSubImage(IplImage **img)
{
	/* release img2 when you don't need it any more.*/
	cvReleaseImage(img);
}

/*
save roi to a file
*/
FILE *store_roi(int reset, IplImage *img, CvRect roi)
{
	static int no;
	static struct timeval t0;
	static char subdir[80];
	static uint64_t mt0=0;
	static FILE *flist=NULL;
	struct timeval t1;
	uint64_t mt1;
	unsigned int e;
	char fstr[256];
	int i,j;
	float rs=0.0,r=0.0,g=0.0,b=0.0;

	printf("store_roi (%d,%d,%d,%d)\n",roi.x,roi.y,roi.width,roi.height);

	if(reset) no=0;
	if(!no){
		//init
		struct stat s;
		struct timeval  tv;
	   	struct tm *tm0;
		gettimeofday(&tv, NULL);
    	if((tm0 = localtime(&tv.tv_sec)) == NULL){
    		perror("localtime error!\n");
    		return NULL;
    	}
   		strftime (subdir, 80, "%m%d%H%M%S", tm0);
		printf("subdir:[%s] %s\n", subdir, asctime(tm0));
		if (stat(subdir, &s) == -1) {
			if(mkdir(subdir, 0722)){
				printf("folder %s create failed.\n", subdir);
				return NULL;
			}
		}
		gettimeofday(&t0, NULL);
		mt0 = (t0.tv_sec * 1000) + t0.tv_usec/1000;
		//create a file list.txt
		if(flist) fclose(flist);
		sprintf(fstr, "%s/list.txt",subdir);
		flist=fopen(fstr,"w");
		if(flist == NULL){
			printf("%s error\n",fstr);
			return NULL;
		}
		printf("The IplImage Size is : %d\n",sizeof(IplImage));
		printf("The nSize            is : %d\n",img->nSize);
		printf("The width           is : %d\n",img->width);
		printf("The height          is : %d\n",img->height);
		printf("The nChannels     is : %d\n",img->nChannels);
		printf("The depth           is : %d\n",img->depth);
		printf("The widthStep     is : %d\n",img->widthStep);
		printf("The imageSize     is : %d\n",img->imageSize);
		printf("The dataOrder     is : %d\n",img->dataOrder);
		printf("The origin           is : %d\n",img->origin);
		printf("The align            is : %d\n",img->align);
	}
	gettimeofday(&t1, NULL);
	mt1 = (t1.tv_sec * 1000) + t1.tv_usec/1000;
	e = mt1 - mt0;
	//printf("mt0=%lu, mt1=%lu, e=%u\n", mt0, mt1, e);
	//file name : f%04u_ticks.bmp
	sprintf(fstr, "%s/f%04u-%08u.bmp",subdir,no,e);

	printf("storing:%s\n",fstr);
	/* sets the Region of Interest
	   Note that the rectangle area has to be __INSIDE__ the image */
	cvSetImageROI(img, roi);
	cvSaveImage(fstr,img, 0);
	/* always reset the Region of Interest */
	cvResetImageROI(img);

	/* averaging ROI's pixel*/
	for(i=roi.y;i< (roi.y+roi.height);i++)
    {
        for(j=roi.x*img->nChannels;j< (roi.x + roi.width) * img->nChannels;j+=img->nChannels)
        {
            b += (unsigned char)img->imageData[i*img->widthStep+j];
            g += (unsigned char)img->imageData[i*img->widthStep+j+1];
            r += (unsigned char)img->imageData[i*img->widthStep+j+2];
        }
    }
    rs=(roi.width*roi.height);
    //update the bmp file name to the list.txt
	//printf("rs=%.1f, bgr=%.1f %.1f,%.1f %.1f,%.1f %.1f\n",rs, b,b/rs, g,g/rs, r,r/rs);
	fprintf(flist, "f%04u-%08u.bmp %.3f %.3f %.3f\n",no,e,b/rs, g/rs, r/rs);

	no++;

	return flist;
}

void process_image(IplImage* frame, int draw, int print)
{
	int i;

	unsigned char* pixels;
	int nrows, ncols, ldim;

	int ndetections;
	float qs[MAXNDETECTIONS], rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS];

	static struct timeval t1;
	struct timeval t2;
    unsigned int e;

	static IplImage* gray = 0;

	// a structure that encodes object appearance
	static char appfinder[] =
		{
			#include "../../cascades/facefinder.ea"
		};

	// grayscale image
	if(!gray)
		gray = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, 1);
	if(frame->nChannels == 3)
		cvCvtColor(frame, gray, CV_RGB2GRAY);
	else
		cvCopy(frame, gray, 0);

	// get relevant image data
	pixels = (unsigned char*)gray->imageData;
	nrows = gray->height;
	ncols = gray->width;
	ldim = gray->widthStep;

	// actually, all the smart stuff happens here
#ifndef _ROTATION_INVARIANT_DETECTION_
	ndetections = find_objects(0.0f, rs, cs, ss, qs, MAXNDETECTIONS, appfinder, pixels, nrows, ncols, ldim, SCALEFACTOR, STRIDEFACTOR, minsize, MIN(nrows, ncols), 1);
#else
	// scan the image at 12 different orientations
	ndetections = 0;
	#pragma omp parallel for
	{
		for(i=0; i<12; ++i)
		{
			float orientation = i*2*3.14f/12;// 2Pi/12= 360"/12=30"

			ndetections += find_objects(orientation, &rs[ndetections], &cs[ndetections], &ss[ndetections], &qs[ndetections], MAXNDETECTIONS-ndetections, appfinder, pixels, nrows, ncols, ldim, SCALEFACTOR, STRIDEFACTOR, minsize, MIN(nrows, ncols), 1);
		}
	}
#endif

	// if the flag is set, draw each detection
	if(draw){
		//frame per ms
		gettimeofday(&t2, NULL);
		uint64_t mt2 = (t2.tv_sec * 1000) + t2.tv_usec/1000;
		uint64_t mt1 = (t1.tv_sec * 1000) + t1.tv_usec/1000;
		e = mt2 - mt1;
		//printf("\nt1=(%d-%d), t2=(%d-%d) , e=%lums\n", t1.tv_sec, t1.tv_usec, t2.tv_sec, t2.tv_usec, e);
		printf("fps=%3.1f e=%3ums,\n", 1000.0/e, e);
		t1 = t2;

		for(i=0; i<ndetections; ++i){
			if(qs[i]>=QCUTOFF){ // check the confidence threshold
				//void cvCircle(CvArr* img, CvPoint center, int radius, CvScalar color, int thickness=1, int line_type=8, int shift=0 )
				//cvCircle(frame, cvPoint(cs[i], rs[i]), ss[i]/2, CV_RGB(255, 0, 0), 4, 8, 0); // we draw circles here since height-to-width ratio of the detected face regions is 1.0f

				/*orginal face ROI : (x-D/2, y-D/2) - (x+D/2, y+D/2)
				  get rid of top hair area of the face : D/5
				  The ROI to send depedning on the ROI area. If the area changes, the ratio should be changed to
				  cover the most of skin area.
				*/
				int xoff,yoff1, yoff2, d;
				d = (int)ss[i];
				xoff = d*0.3;
				yoff2 = yoff1 = d>>1;
				yoff1 -= d/4;
				//hair
				cvRectangle(frame, cvPoint(cs[i]- xoff, rs[i]- (d>>1)), cvPoint(cs[i]+ xoff, rs[i]-yoff1), CV_RGB(255, 0, 0), 4, 8, 0);
				//face raw trace region which is smaller than face ROI.
				cvRectangle(frame, cvPoint(cs[i]-xoff, rs[i]-yoff1), cvPoint(cs[i]+xoff, rs[i]+yoff2), CV_RGB(0, 255, 0), 4, 8, 0);
			}
		}
	}
	// if the flag is set, print the results to standard output
	if(print)
	{
		for(i=0; i<ndetections; ++i)
			if(qs[i]>=QCUTOFF) // check the confidence threshold
				printf("%d %d %d %f\n", (int)rs[i], (int)cs[i], (int)ss[i], qs[i]);
	}
}

void process_webcam_frames(int idx)
{
	static uint64_t ut1;
	CvCapture* capture;
	struct v4l2_frmivalenum frmival;
	uint32_t fps;
	IplImage* frame=NULL;
	IplImage* framecopy=NULL;
	int save_roi=0;
	int stop;
	int camfd;
	struct timeval pt2;
	uint64_t ut2;
	const char* windowname = "webcam";

	//cvCaptureFromCAM will lock camera device, so do VIDIOC_S_FMT before openCV.
	preset_videofmt(idx);
	// try to initialize video capture from the default webcam
	capture = cvCaptureFromCAM(idx);
	if(!capture)
	{
		printf("Cannot initialize video capture!\n");
		return;
	}
	camfd = cvGetCamFD(capture);
	postset_cam(camfd);

	// start the main loop in which we'll process webcam output
	framecopy = 0;
	stop = 0;
	cvNamedWindow(windowname,CV_WINDOW_AUTOSIZE);
	while(!stop)
	{
		static FILE *file=NULL;

		// retrieve a pointer to the image acquired from the webcam
		if(!cvGrabFrame(capture))
			break;
		frame = cvRetrieveFrame(capture, 1);
		//processing is done
		gettimeofday(&pt2, NULL);
		ut2 = (pt2.tv_sec * 1000000) + pt2.tv_usec;
		if( ut1 && (ut2 > ut1)){
//			printf("\npt=%lu us, fps=%.1f\n", ut2-ut1, 1000000.0/(ut2-ut1));
			printf("fps=%.0f\n", 1000000.0/(ut2-ut1));
		}
		ut1=ut2;

		// wait 5 miliseconds
		int key = cvWaitKey(1);
		if(key == 's'){
			save_roi=1;
			printf("save_roi=1\n");
		}else if(key == 'e'){
			save_roi=0;
			if(file)
				fclose(file);
			file=NULL;
			printf("save_roi=0\n");
		}

		// we terminate the loop if we don't get any data from the webcam or the user has pressed 'q'
		if(!frame || key=='q')
			stop = 1;
		else
		{
			CvRect roi,sroi;
			// we mustn't tamper with internal OpenCV buffers and that's the reason why we're making a copy of the current frame
			//framecopy = cvCloneImage(frame);
			if(!framecopy)
				framecopy = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, frame->nChannels);
			cvCopy(frame, framecopy, 0);

			// webcam outputs mirrored frames (at least on my machines); you can safely comment out this line if you find it unnecessary
			//cvFlip(framecopy, framecopy, 1);

			// all the smart stuff happens in the following function
			//process_image(framecopy, 1, 0);

			//predefined ROI (2w/5,2h/5) - (2.5W/4, 2.5h/4)
			roi.x = (frame->width - ROI_BORDERW*2 - ROI_WIDTH)/2;
			roi.y = (frame->height - ROI_BORDERW*2 - ROI_HEIGHT)/2;
			roi.width = ROI_WIDTH + ROI_BORDERW*2;
			roi.height = ROI_HEIGHT + ROI_BORDERW*2 ;
			/*
			roi.x = frame->width*1.8/5.0;
			roi.y = frame->height*1.0/5.0;
			roi.width = frame->width*1.4/5.0;
			roi.height = frame->height*2.25/5.0;
			*/
			sroi=roi;
			sroi.x += (ROI_BORDERW);
			sroi.y += (ROI_BORDERW);
			sroi.width = ROI_WIDTH;
			sroi.height = ROI_HEIGHT;
			//printf(" roi=%d %d %d %d\n", roi.x, roi.y, roi.width,roi.height);
			//printf("sroi=%d %d %d %d\n", sroi.x, sroi.y, sroi.width,sroi.height);
			cvRectangle(framecopy, cvPoint(roi.x, roi.y),
						cvPoint(roi.x+roi.width, roi.y+roi.height), CV_RGB(255, 0, 0), ROI_BORDERW,8, 0);
			//printf("save_roi=[%d]\n",save_roi );
			if(save_roi)
				file=store_roi(0, framecopy, sroi);
			// display the image to the user
			cvShowImage(windowname, framecopy);
		}
	}
exit:
	// cleanup
	if(framecopy)
		cvReleaseImage(&framecopy);
	if(frame)
		cvReleaseImage(&frame);

	cvReleaseCapture(&capture);
	cvDestroyWindow(windowname);
}

int main(int argc, char* argv[])
{
	int fd, i=0;
	char dev[20];
	if(argc == 2)
		sscanf(argv[1], "%d", &i);//camera index from 0
	process_webcam_frames(i);
#if 0
	preset_videofmt(i);
	sprintf(dev,"/dev/video%d",i);
	printf("cam%d:%s\n",i, dev);
    fd = open(dev, O_RDWR);
    if (fd == -1){
		perror("Opening video device");
        return 1;
    }
	postset_cam(fd);
	init_mmap(fd);
	while(1){
		int key = cvWaitKey(1);
		capture_image(fd,"windows");
		// wait 5 miliseconds
		if(key == 'q'){
			break;
		}
	}
#endif
	close(fd);

	return 0;
}

