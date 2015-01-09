/*
 *  V4L2 video capture example
 *
 *  This program can be used and distributed without restrictions.
 *
 *	This program were got from V4L2 API, Draft 0.20
 *		available at: http://v4l2spec.bytesex.org/
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>             /* getopt_long() */
#include <fcntl.h>              /* low-level i/o */
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/ioctl.h>

#include <linux/videodev2.h>

#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "dsp-ica.h"

#define FORCED_WIDTH  640
#define FORCED_HEIGHT 480
#define FORCED_FORMAT V4L2_PIX_FMT_YUYV	//V4L2_PIX_FMT_MJPEG
#define FORCED_FIELD  V4L2_FIELD_ANY
#define	ROI_BORDERW		(4)
#define	ROI_WIDTH		(160)
#define	ROI_HEIGHT		(180)
#define ROI_Y_OFFSET	(90)

/*http://forum.processing.org/one/topic/webcam-with-stable-framerate-25fps-on-high-resolution.html
The stable fps can only occur in the low frame rate. Higher fps>15 may not be stable.
*/
//#define FORCED_FPS		(9)	//fps9 or fps10 are the most stable fps.

#define CLEAR(x) memset(&(x), 0, sizeof(x))

enum io_method {
	IO_METHOD_READ,
	IO_METHOD_MMAP,
	IO_METHOD_USERPTR,
};

struct buffer {
	void   *start;
	size_t  length;
};

static PCircularQueue cq,cqLong;

static char            *dev_name;
static enum io_method   io = IO_METHOD_MMAP;
static int              fd = -1;
struct buffer          *buffers;
static unsigned int     n_buffers;
static int				out_buf;
static int              force_format;
static int              frame_count = 0;

static char *windowname="v4l2 capture";

static void deinitCQs(void)
{
	pr_debug(DSP_INFO,"+%s:\n", __func__);

	if(cqLong) cqDeinit(cqLong);
	cqLong=NULL;
	if(cq) cqDeinit(cq);
	cq=NULL;
}

static void initCQs(void)
{
	pr_debug(DSP_INFO,"+%s:\n", __func__);
	deinitCQs();
	cq=cqInit(getSampleWindow() * getFPS());
	cqLong=cqInit(getMaxSampleTime() * getFPS());
}

static void errno_exit(const char *s)
{
	pr_debug(DSP_ERROR, "+%s:%s error %d, %s\n", __func__, s, errno, strerror(errno));
	deinitCQs();
	exit(EXIT_FAILURE);
}

static int xioctl(int fh, int request, void *arg)
{
	int r;

	do {
		r = ioctl(fh, request, arg);
	} while (-1 == r && EINTR == errno);

	return r;
}

int QueryCap(int fd, struct v4l2_capability *pcaps)
{
    struct v4l2_capability caps = {0};
    if (-1 == xioctl(fd, VIDIOC_QUERYCAP, pcaps)){
        perror("Querying Capabilities");
        return 1;
    }
    return 0;
}

int QueryCropCap(int fd, struct v4l2_cropcap *pcropcap)
{
    pcropcap->type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    if (-1 == xioctl (fd, VIDIOC_CROPCAP, pcropcap)) {
        perror("Querying Cropping Capabilities");
        return 1;
    }
}

int GetAutoWhiteBalance(int fd)
{
	struct v4l2_control ctrl ={0};
	ctrl.id = V4L2_CID_AUTO_WHITE_BALANCE;
	if (-1 == xioctl(fd, VIDIOC_G_CTRL, &ctrl)){
		perror("getting V4L2_CID_AUTO_WHITE_BALANCE");
	}
//	printf("V4L2_CID_AUTO_WHITE_BALANCE = 0x%x\n",ctrl.value );
	pr_debug(DSP_DEBUG,"-%s:%d\n", __func__, ctrl.value);
	return ctrl.value;
}

int SetAutoWhiteBalance(int fd, int enable)
{
	struct v4l2_control ctrl ={0};
	pr_debug(DSP_DEBUG,"%s:\n", __func__);
	ctrl.id = V4L2_CID_AUTO_WHITE_BALANCE;
	ctrl.value = enable;
	if (-1 == xioctl(fd, VIDIOC_S_CTRL, &ctrl)){
		perror("setting V4L2_CID_AUTO_WHITE_BALANCE");
	}
//	printf("V4L2_CID_AUTO_WHITE_BALANCE = (0x%x -> 0x%x)\n",enable, GetAutoWhiteBalance(fd) );
	return GetAutoWhiteBalance(fd) == enable;
}

/* v4l2-controls.h
enum  v4l2_exposure_auto_type {
	V4L2_EXPOSURE_AUTO = 0,
	V4L2_EXPOSURE_MANUAL = 1,
	V4L2_EXPOSURE_SHUTTER_PRIORITY = 2,
	V4L2_EXPOSURE_APERTURE_PRIORITY = 3
};
V4L2_CID_EXPOSURE_AUTO 	enum v4l2_exposure_auto_type
 	Enables automatic adjustments of the exposure time and/or iris aperture.
 	The effect of manual changes of the exposure time or iris aperture while
 	these features are enabled is undefined, drivers should ignore such requests.
 	Possible values are:
	V4L2_EXPOSURE_AUTO 	Automatic exposure time, automatic iris aperture.
	V4L2_EXPOSURE_MANUAL 	Manual exposure time, manual iris.
	V4L2_EXPOSURE_SHUTTER_PRIORITY 	Manual exposure time, auto iris.
	V4L2_EXPOSURE_APERTURE_PRIORITY 	Auto exposure time, manual iris.

V4L2_CID_EXPOSURE_ABSOLUTE 	integer
 	Determines the exposure time of the camera sensor. The exposure time is
 	limited by the frame interval. Drivers should interpret the values as 100 µs
 	units, where the value 1 stands for 1/10000th of a second, 10000 for
 	1 second and 100000 for 10 seconds.

V4L2_CID_EXPOSURE_AUTO_PRIORITY 	boolean
 	When V4L2_CID_EXPOSURE_AUTO is set to AUTO or APERTURE_PRIORITY, this
 	control determines if the device may dynamically vary the frame rate.
 	By default this feature is disabled (0) and the frame rate must remain constant.

V4L2_CID_EXPOSURE_BIAS 	integer menu
 	Determines the automatic exposure compensation, it is effective only
 	when V4L2_CID_EXPOSURE_AUTO control is set to AUTO, SHUTTER_PRIORITY or
 	APERTURE_PRIORITY. It is expressed in terms of EV, drivers should interpret
 	the values as 0.001 EV units, where the value 1000 stands for +1 EV.
	Increasing the exposure compensation value is equivalent to decreasing the
	exposure value (EV) and will increase the amount of light at the image
	sensor. The camera performs the exposure compensation by adjusting absolute
	exposure time and/or aperture.

V4L2_CID_EXPOSURE_METERING 	enum v4l2_exposure_metering
 	Determines how the camera measures the amount of light available for the
 	frame exposure. Possible values are:

	V4L2_EXPOSURE_METERING_AVERAGE 	Use the light information coming from the
		entire frame and average giving no weighting to any particular portion of
		the metered area.
	V4L2_EXPOSURE_METERING_CENTER_WEIGHTED 	Average the light information coming
		from the entire frame giving priority to the center of the metered area.
	V4L2_EXPOSURE_METERING_SPOT 	Measure only very small area at the center
		of the frame.
	V4L2_EXPOSURE_METERING_MATRIX 	A multi-zone metering. The light intensity
		is measured in several points of the frame and the the results are
		combined. The algorithm of the zones selection and their significance
		in calculating the final value is device dependent.
*/
static void autoExposureType(int type)
{
   	switch(type){
   	case V4L2_EXPOSURE_AUTO:
   		pr_debug(DSP_INFO,"V4L2_EXPOSURE_AUTO[%d]\n",type);
   		break;
   	case V4L2_EXPOSURE_MANUAL:
   		pr_debug(DSP_INFO,"V4L2_EXPOSURE_MANUAL[%d]\n",type);
   		break;
   	case V4L2_EXPOSURE_SHUTTER_PRIORITY:
   		pr_debug(DSP_INFO,"V4L2_EXPOSURE_SHUTTER_PRIORITY[%d]\n",type);
   		break;
   	case V4L2_EXPOSURE_APERTURE_PRIORITY:
   		pr_debug(DSP_INFO,"V4L2_EXPOSURE_APERTURE_PRIORITY[%d]\n",type);
   		break;
   	default:
   		perror("unknown exposure type");
   	}
}

/*
V4L2_EXPOSURE_MANUAL and V4L2_EXPOSURE_APERTURE_PRIORITY are commonly used.
*/
int SetAutoExposure(int fd, int type)
{
	struct v4l2_control ctrl ={0};
	ctrl.id = V4L2_CID_EXPOSURE_AUTO;
   	ctrl.value = type;
   	pr_debug(DSP_INFO,"SetAutoExposure=");
   	autoExposureType(type);
   	if (-1 == xioctl(fd,VIDIOC_S_CTRL,&ctrl)) {
      perror("setting V4L2_CID_EXPOSURE_AUTO");
      return -1;
	}
	return type;
}

int GetAutoExposure(int fd)
{
	struct v4l2_control ctrl ={0};
	ctrl.id = V4L2_CID_EXPOSURE_AUTO;
	if (-1 == xioctl(fd, VIDIOC_G_CTRL, &ctrl)){
		perror("getting V4L2_CID_EXPOSURE_AUTO");
		return -1;
	}
	pr_debug(DSP_INFO,"GetAutoExposure=");
	autoExposureType(ctrl.value);
	return ctrl.value;
}

int SetAutoExposureAutoPriority(int fd, int p)
{
	struct v4l2_control ctrl ={0};
	ctrl.id = V4L2_CID_EXPOSURE_AUTO_PRIORITY ;
   	ctrl.value = p;
   	if (-1 == xioctl(fd,VIDIOC_S_CTRL,&ctrl)) {
      perror("setting V4L2_CID_EXPOSURE_AUTO_PRIORITY");
      return -1;
	}
	return p;
}

int GetAutoExposureAutoPriority(int fd)
{
	struct v4l2_control ctrl ={0};
	ctrl.id = V4L2_CID_EXPOSURE_AUTO_PRIORITY ;
   	if (-1 == xioctl(fd,VIDIOC_G_CTRL,&ctrl)) {
      perror("getting V4L2_CID_EXPOSURE_AUTO_PRIORITY");
      return -1;
	}
	return ctrl.value;
}

/*
	The exposure time is limited by the frame interval. Drivers should interpret
	the values as 100 µs units, where the value 1 stands for 1/10000th of a second,
	10000 for 1 second and 100000 for 10 seconds.
	if 30fps is the goal, each frame takes 1/30s = 33ms = 33*10*100us =
	330 * 100us. It depends on the ambient light to tuning exposure time.
*/
int SetManualExposure(int fd, int val)
{
	struct v4l2_control ctrl ={0};
	if(SetAutoExposure(fd, V4L2_EXPOSURE_MANUAL) != V4L2_EXPOSURE_MANUAL){
		perror("setting V4L2_EXPOSURE_MANUAL");
		return -1;
	}

	ctrl.id = V4L2_CID_EXPOSURE_ABSOLUTE;
   	ctrl.value = val;
   	if (-1 == xioctl(fd,VIDIOC_S_CTRL,&ctrl)) {
      perror("setting V4L2_CID_EXPOSURE_ABSOLUTE");
      return -1;
	}
	return val;
}

int GetManualExposure(int fd)
{
	struct v4l2_control ctrl ={0};
/*	if(GetAutoExposure(fd) != V4L2_EXPOSURE_MANUAL){
		perror("SetManualExposure is not in V4L2_EXPOSURE_MANUAL");
		return -1;
	}
*/
	ctrl.id = V4L2_CID_EXPOSURE_ABSOLUTE;
   	if (-1 == xioctl(fd,VIDIOC_G_CTRL,&ctrl)) {
      perror("getting V4L2_EXPOSURE_MANUAL");
      return -1;
	}
	printf("GetManualExposure=%d\n",ctrl.value);
	return ctrl.value;
}

int EnumVideoFMT(int fd)
{
	int support_grbg10 = 0;
	struct v4l2_fmtdesc fmtdesc = {0};
    fmtdesc.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    char fourcc[5] = {0};
    char c, e;
    printf("\n  FMT : CE Desc\n--------------------\n");
    while (0 == xioctl(fd, VIDIOC_ENUM_FMT, &fmtdesc))
    {
        strncpy(fourcc, (char *)&fmtdesc.pixelformat, 4);
        if (fmtdesc.pixelformat == V4L2_PIX_FMT_SGRBG10)
            support_grbg10 = 1;
        c = fmtdesc.flags & 1? 'C' : ' ';
        e = fmtdesc.flags & 2? 'E' : ' ';
        printf("  %s: [%c][%c], [%s]\n", fourcc, c, e, fmtdesc.description);
        fmtdesc.index++;
    }
}

int GetVideoFMT(int fd, struct v4l2_format *pfmt)
{
    int i;
//    struct v4l2_format fmt = {0};
    char fourcc[5] = {0};
    pfmt->type = V4L2_BUF_TYPE_VIDEO_CAPTURE; //V4L2_BUF_TYPE_VIDEO_CAPTURE_MPLANE;
	for(i = 0; i < 20; i ++){
		if (-1 == xioctl(fd, VIDIOC_G_FMT, pfmt))
		{
		    perror("getting Pixel Format");
			usleep(5000);
			continue;
		}else{
		strncpy(fourcc, (char *)&pfmt->fmt.pix.pixelformat, 4);
		printf( "Gotten Camera Mode:\n"
		    "  Width: %d\n"
		    "  Height: %d\n"
		    "  PixFmt: %s\n"
		    "  Field: %d\n"
		    "  priv: 0x%x\n",
		    pfmt->fmt.pix.width,
		    pfmt->fmt.pix.height,
		    fourcc,
		    pfmt->fmt.pix.field,
		    pfmt->fmt.pix.priv);
		    break;
		}
	}
}

/*
http://permalink.gmane.org/gmane.linux.drivers.uvc.devel/4537
You can't change the format, VIDIOC_S_FMT, while buffers are allocated.
You need to free the buffers before, using VIDIOC_REQBUFS with a buffer count of 0.
*/
int SetVideoFMT(int fd, struct v4l2_format fmt)
{
    int i;
//    struct v4l2_format fmt = {0};
    char fourcc[5] = {0};
    fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
//    fmt.fmt.pix.width = 640;
//    fmt.fmt.pix.height = 480;
    //fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_BGR24;
    //fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_GREY;
//    fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_MJPEG;
//    fmt.fmt.pix.pixelformat = pixelformat;//V4L2_PIX_FMT_YUYV;
//    fmt.fmt.pix.field = V4L2_FIELD_NONE;
//	for(i = 0; i < 20; i ++){
		if (-1 == xioctl(fd, VIDIOC_S_FMT, &fmt))
		{
		    perror("Setting Pixel Format");
			//usleep(10000);
			//continue;
		}else{
		strncpy(fourcc, (char *)&fmt.fmt.pix.pixelformat, 4);
		printf( "Selected Camera Mode:\n"
		    "  Width: %d\n"
		    "  Height: %d\n"
		    "  PixFmt: %s\n"
		    "  Field: %d\n",
		    fmt.fmt.pix.width,
		    fmt.fmt.pix.height,
		    fourcc,
		    fmt.fmt.pix.field);
		    //break;
		}
	//}
	return 0;
}

int SetFPSParam(int fd, uint32_t fps)
{
	double fps_new=(double)fps;
	struct v4l2_streamparm param;
    memset(&param, 0, sizeof(param));
    param.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    param.parm.capture.timeperframe.numerator = 1;
    param.parm.capture.timeperframe.denominator = fps;
	if (-1 == xioctl(fd, VIDIOC_S_PARM, &param)){
		perror("unable to change device parameters");
		return ;
	}

	if(param.parm.capture.timeperframe.numerator){
		fps_new = param.parm.capture.timeperframe.denominator
	                 / param.parm.capture.timeperframe.numerator;
		if ((double)fps != fps_new) {
			printf("unsupported frame rate [%d,%f]\n", fps, fps_new);
			return;
		}else{
			printf("new fps:%u , %u/%u\n",fps, param.parm.capture.timeperframe.denominator,
			param.parm.capture.timeperframe.numerator);
		}
	}
}

uint32_t GetFPSParam(int fd, double fps, struct v4l2_frmivalenum *pfrmival)
{
    struct v4l2_frmivalenum frmival[10];
    float fpss[10];
    int i=0;
/*
    memset(&frmival,0,sizeof(frmival));
    frmival.pixel_format = fmt;
    frmival.width = width;
    frmival.height = height;*/
    memset(fpss,0,sizeof(fpss));
    while(pfrmival->index < 10){
	    if (-1 == xioctl(fd, VIDIOC_ENUM_FRAMEINTERVALS, pfrmival)){
		    perror("getting VIDIOC_ENUM_FRAMEINTERVALS");
		    break;
		}
		frmival[pfrmival->index] = *pfrmival;
        if (pfrmival->type == V4L2_FRMIVAL_TYPE_DISCRETE){
	    	double f;
        	f = (double)pfrmival->discrete.denominator/pfrmival->discrete.numerator;
        	pr_debug(DSP_INFO, "[%u/%u]\n", pfrmival->discrete.denominator,
        						pfrmival->discrete.numerator);
            pr_debug(DSP_INFO,"[%dx%d] %f fps\n", pfrmival->width, pfrmival->height,f);

			fpss[pfrmival->index]=f;
			frmival[pfrmival->index]=*pfrmival;
        }else{
        	double f1,f2;
        	f1 = (double)pfrmival->stepwise.max.denominator/pfrmival->stepwise.max.numerator;
        	f2 = (double)pfrmival->stepwise.min.denominator/pfrmival->stepwise.min.numerator;
            pr_debug(DSP_INFO,"[%dx%d] [%f,%f] fps\n", pfrmival->width, pfrmival->height,f1,f2);
       	}
       	pr_debug(DSP_INFO, "idx=%d\n", pfrmival->index);
       	pfrmival->index++;
    }
    /* list is in increasing order */
    if(pfrmival->index){
    	i = pfrmival->index;
	    while(--i >= 0){
    		if(fps <= fpss[i] ){
    			break;
    		}
    	}
    	*pfrmival = frmival[i];
    	pr_debug(DSP_DEBUG, "fps:request[%.1f],but set to [%.1f]\n", fps, fpss[i]);
    }
    return (uint32_t)fpss[i];
}

int PrintFrameInterval(int fd, unsigned int fmt, unsigned int width, unsigned int height)
{
    struct v4l2_frmivalenum frmival;
    memset(&frmival,0,sizeof(frmival));
    frmival.pixel_format = fmt;
    frmival.width = width;
    frmival.height = height;
    while(1){
	    if (-1 == xioctl(fd, VIDIOC_ENUM_FRAMEINTERVALS, &frmival)){
		    perror("getting VIDIOC_ENUM_FRAMEINTERVALS");
		    return -1;
		}

        if (frmival.type == V4L2_FRMIVAL_TYPE_DISCRETE){
           	printf("[%u/%u]\n", frmival.discrete.denominator,
        						frmival.discrete.numerator);
            printf("[%dx%d] %f fps\n", width, height,
            1.0*frmival.discrete.denominator/frmival.discrete.numerator);
		}else
            printf("[%dx%d] [%f,%f] fps\n", width, height,
            1.0*frmival.stepwise.max.denominator/frmival.stepwise.max.numerator,
            1.0*frmival.stepwise.min.denominator/frmival.stepwise.min.numerator);
        frmival.index++;
    }
    return 0;
}

int EnumFrameRate(int fd, unsigned int format)
{
    unsigned int width=0, height=0;;
    struct v4l2_frmsizeenum frmsize;
    memset(&frmsize,0,sizeof(frmsize));
    frmsize.pixel_format = format; //V4L2_PIX_FMT_JPEG;
	printf("\nEnumFrameRate:\n");
    while(1){
	    if (-1 == xioctl(fd, VIDIOC_ENUM_FRAMESIZES, &frmsize)){
		    perror("getting VIDIOC_ENUM_FRAMESIZES");
		    return -1;
		}
		printf("frmsize.type=%d\n", frmsize.type);
        if (frmsize.type == V4L2_FRMSIZE_TYPE_DISCRETE){
            PrintFrameInterval(fd, frmsize.pixel_format, frmsize.discrete.width,
            frmsize.discrete.height);
        }else{
            for (width=frmsize.stepwise.min_width; width< frmsize.stepwise.max_width;
            	width+=frmsize.stepwise.step_width)
                for (height=frmsize.stepwise.min_height;
                	height< frmsize.stepwise.max_height;
                	height+=frmsize.stepwise.step_height)
                    PrintFrameInterval(fd, frmsize.pixel_format, width, height);
        }
        frmsize.index++;
    }
    return 0;
}

int print_caps(int fd)
{
    struct v4l2_capability caps = {0};
    struct v4l2_cropcap cropcap = {0};

    if (!QueryCap(fd, &caps)){
	    printf( "Driver Caps:\n"
            "  Driver: \"%s\"\n"
            "  Card: \"%s\"\n"
            "  Bus: \"%s\"\n"
            "  Version: %d.%d\n"
            "  Capabilities: %08x\n",
            caps.driver,
            caps.card,                caps.bus_info,
            (caps.version>>16)&&0xff,
            (caps.version>>24)&&0xff,
            caps.capabilities);
	}

    if (!QueryCropCap(fd, &cropcap)){
    	printf( "Camera Cropping:\n"
        "  Bounds: %dx%d+%d+%d\n"
        "  Default: %dx%d+%d+%d\n"
        "  Aspect: %d/%d\n",
        cropcap.bounds.width, cropcap.bounds.height, cropcap.bounds.left, cropcap.bounds.top,
        cropcap.defrect.width, cropcap.defrect.height, cropcap.defrect.left, cropcap.defrect.top,
        cropcap.pixelaspect.numerator, cropcap.pixelaspect.denominator);
	}
    EnumVideoFMT(fd);
    //int support_grbg10 = 0;
    /*
    if (!support_grbg10)
    {
        printf("Doesn't support GRBG10.\n");
        return 1;
    }*/

	GetAutoWhiteBalance(fd);

    return 0;
}

/* openCV cvCaptureFromCAM will lock the /dev/videoxxx, so the VIDIOC_S_FMT
can't work.
Set the capture format before openCV open the camera device.
*/
int extra_cam_setting(int camfd)
{
	struct v4l2_format fmt;
	struct v4l2_frmivalenum frmival;
	uint32_t fps;

	print_caps(camfd);
	GetVideoFMT(camfd, &fmt);
	EnumFrameRate(camfd, V4L2_PIX_FMT_YUYV);

	memset(&frmival,0,sizeof(frmival));
    frmival.pixel_format = V4L2_PIX_FMT_YUYV;
    frmival.width = FORCED_WIDTH;
    frmival.height = FORCED_HEIGHT;
	//fps = GetFPSParam(camfd, (double)FORCED_FPS, &frmival);
	fps = GetFPSParam(camfd, (double)getFPS(), &frmival);
	setFPS(fps);	//set the proper fps by the camera's cap
    frmival.pixel_format = V4L2_PIX_FMT_YUYV;
    frmival.width = FORCED_WIDTH;
    frmival.height = FORCED_HEIGHT;
	SetFPSParam(camfd, fps);
	SetAutoWhiteBalance(camfd, 0);
	SetAutoExposure(camfd, V4L2_EXPOSURE_AUTO);
	SetAutoExposureAutoPriority(camfd,0);
#if 0
	SetAutoExposure(camfd, /*V4L2_EXPOSURE_MANUAL ,*/ V4L2_EXPOSURE_APERTURE_PRIORITY  );
	SetAutoExposureAutoPriority(camfd,0);
	printf("AutoPriority=%d\n",GetAutoExposureAutoPriority(camfd));
//	PrintFrameInterval(camfd, frmival.pixel_format, frmival.width, frmival.height);
//	SetAutoExposure(camfd, V4L2_EXPOSURE_MANUAL);
	SetManualExposure(camfd, 110);
	GetManualExposure(camfd);
#endif
}

void cvPrintf(IplImage* img, char *text, CvPoint TextPosition, CvFont Font1,
			  CvScalar Color)
{
  //char text[] = "Funny text inside the box";
  //int fontFace = FONT_HERSHEY_SCRIPT_SIMPLEX;
  //double fontScale = 2;

  // center the text
  //Point textOrg((img.cols - textSize.width)/2,
	//			(img.rows + textSize.height)/2);

  //double Scale=2.0;
  //int Thickness=2;
  //CvScalar Color=CV_RGB(255,0,0);
  //CvPoint TextPosition=cvPoint(400,50);
  //CvFont Font1=cvFont(Scale,Thickness);
  // then put the text itself
  cvPutText(img, text, TextPosition, &Font1, Color);
  //cvPutText(img, text, TextPosition, &Font1, CV_RGB(0,255,0));
}

int spectraAnalysis(IplImage* framecopy, CvScalar rgb, float *pr, float *rr)
{
	int i, n=0;
	pr_debug(DSP_INFO,"+%s:\n", __func__);
	//PCircularQueue pq[]={cq, cqLong};
	PCircularQueue pq[]={cq};
	int win_samples[]={getSampleWindow()*getFPS(), getMaxSampleTime() * getFPS()};
	//insert to cq and process it
	for(i = 0; i < sizeof(pq)/sizeof(PCircularQueue); i ++){
		if(cqInsert(pq[i], (PElemType)&rgb)){
			pr_debug(DSP_DEBUG,"!!!!queue[%d] is full and is wating...\n", i);
		}else{
			int j;
			int steps=getStepping()*getFPS();
			pr_debug(DSP_INFO,"cq[%d]: %d\n",i, cqUsedSize(pq[i]));
			//indicator to count bpm
			char hrtext[20];
			sprintf(hrtext,"[%d/%d]", cqUsedSize(pq[i])+1,cqSize(pq[i]));
			cvPrintf(framecopy, hrtext, cvPoint(200,90), cvFont(2.0,1.0), CV_RGB(0,255,0));

			if(cqFull(pq[i])){
				int row,col;
				CvMat *rawTrace;
				//read queue to rawtrace
				pr_debug(DSP_DEBUG,"###cq[%d] is full : %d, process it\n", i,
						 cqUsedSize(pq[i]));
				pr_debug(DSP_INFO,"N=%d\n", win_samples[i]);
				//3 rows (BGR) * N
				rawTrace = cvCreateMat(MAX_CHANNELS, win_samples[i], RAW_ELEMENT_TYPE);
				if(rawTrace){
					CvScalar spectra[MAX_SPECTRAS];//(magnitude, freq, channel)
					cqClone(pq[i], rawTrace);
					n=rawTraceSpectra(rawTrace, spectra);
					//clear out steps from cq and waits for when cq is full.
					for(j=0;j<steps;j++) cqDelete(pq[i], NULL);
					pr_debug(DSP_INFO,"[%d] frames are cleared out: %d left\n", steps,
							cqUsedSize(pq[i]));
					if(n && pr)
					  *pr=spectra[0].val[1];
					else n=-1;//no valid spectra is found!
					for(j=0; j<n;j++) {
						pr_debug(DSP_DEBUG,"\n***[%d: %.1fbpm (%.3fhz), mag=%.3f,"
								 "ch=%.0f]\n\n",
								 j, spectra[j].val[1]*60.0,spectra[j].val[1],
								spectra[j].val[0],spectra[j].val[2]);
					}
					cvReleaseMat(&rawTrace);
				}
			}
		}
	}
  pr_debug(DSP_INFO,"-%s:n=%d, freq=%.3f\n", __func__,n,pr?*pr:0);
  return n;
}

/*
p is a YUYV 422 format, so 640x480x16bits = 61440 bytes
*/
static CvScalar processFrame(const void *p, int size)
{
	IplImage* framecopy=NULL;
	static uint64_t ut1;
	uint64_t ut2;
	struct timeval pt2;
	pr_debug(DSP_INFO,"%s: size=0x%x\n", __func__, size);

#if 0 	//V4L2_PIX_FMT_MJPEG
    IplImage* frame;
    CvMat cvmat = cvMat(FORCED_HEIGHT, FORCED_WIDTH, CV_8UC3, (void*)p);//MJPEG buffer p
    frame = cvDecodeImage(&cvmat, 1);//sometimes a corrupted frame is retrieved,
    								//so frame should be checked if it's null.
	if(frame && !framecopy){
//		CvSize s = cvGetSize(frame);
//		printf("[%d,%d]\n",s.width, s.height);
		framecopy = cvCreateImage(cvGetSize(frame), IPL_DEPTH_8U, 3);
//		printf("cvCreateImage=%p\n", framecopy);
	}
	if(frame){
		cvCopy(frame, framecopy, 0);
//    	cvShowImage(windowname, frame);
    	cvShowImage("framecopy", framecopy);
    }else{
	    printf("frame NULL, size=%d\n", size);
    }
//	cvReleaseImage(&frame);
#else
//V4L2_PIX_FMT_YUYV
//	if(size < (FORCED_WIDTH*FORCED_HEIGHT*2)){
//		printf("size too small\n");
//		return ;
//	}
	framecopy = cvCreateImage(cvSize(FORCED_WIDTH,FORCED_HEIGHT), IPL_DEPTH_8U, 3);
	yuyv_to_rgb24(FORCED_WIDTH,FORCED_HEIGHT, (unsigned char *)p, framecopy->imageData);
	float m_roi[3]={0.0};
	CvRect sroi;
	sroi.x = (framecopy->width - ROI_WIDTH)/2 - 1;
	sroi.y = (framecopy->height - ROI_HEIGHT)/2 - 1+ ROI_Y_OFFSET;
	sroi.width = ROI_WIDTH;
	sroi.height = ROI_HEIGHT ;
	roi_mean(framecopy, sroi, m_roi);
	pr_debug(DSP_INFO,"m_roi(%.4f,%.4f,%.4f)\n", m_roi[0],m_roi[1],m_roi[2]);
	//
	//draw the bouding ROI on the frame
	/*cvRectangle(framecopy, cvPoint(sroi.x, sroi.y),
				cvPoint(sroi.x+sroi.width, sroi.y+sroi.height),
				CV_RGB(255, 0, 0), ROI_BORDERW,8, 0);*/
#endif
	CvScalar m_rgb=cvScalar(m_roi[0],m_roi[1],m_roi[2], 0.0);
	char hrtext[20];
	float pr=0.0;
	static float valid_pr=0.0;
	int n= spectraAnalysis(framecopy, m_rgb, &pr,NULL);
	static CvScalar textColor={100, 100, 100};//B,G,R
	if(n>0){
	  //green color
	  textColor=CV_RGB(0,255,0);
	  valid_pr=pr;
	}else if(n < 0) {
	  //red color
	  textColor=CV_RGB(255,0,0);
	}
	sprintf(hrtext,"%.1fbpm (%.2fHz)", valid_pr*60, valid_pr);
	cvPrintf(framecopy, hrtext, cvPoint(200,40), cvFont(2.0,2.5), textColor);

	//draw the bouding ROI on the frame
	cvRectangle(framecopy, cvPoint(sroi.x, sroi.y),
				cvPoint(sroi.x+sroi.width, sroi.y+sroi.height),
				textColor, ROI_BORDERW,8, 0);

	gettimeofday(&pt2, NULL);
	ut2 = (pt2.tv_sec * 1000000) + pt2.tv_usec;
	if( ut1 && (ut2 > ut1)){
//			printf("\npt=%lu us, fps=%.1f\n", ut2-ut1, 1000000.0/(ut2-ut1));
	  pr_debug(DSP_INFO,"fps=%.0f\n", 1000000.0/(ut2-ut1));
	  sprintf(hrtext,"fps:%.1f", 1000000.0/(ut2-ut1));
	  cvPrintf(framecopy, hrtext, cvPoint(20,50), cvFont(1.5,1.0), CV_RGB(0,0,255));
	}
	ut1=ut2;

	sprintf(hrtext,"press 'q' to quit...");
	cvPrintf(framecopy, hrtext, cvPoint(10,460),cvFont(1.5,1.0), CV_RGB(255,255, 255));

	cvShowImage(windowname, framecopy);

	if(framecopy)
	  cvReleaseImage(&framecopy);

	return m_rgb;
}

static int readFrame(void)
{
	struct v4l2_buffer buf;
	unsigned int i;
	CvScalar m_rgb=cvScalarAll(0.0f);

	pr_debug(DSP_INFO,"%s:\n", __func__);

	switch (io) {
	case IO_METHOD_READ:
		if (-1 == read(fd, buffers[0].start, buffers[0].length)) {
			switch (errno) {
			case EAGAIN:
				return 0;

			case EIO:
				/* Could ignore EIO, see spec. */

				/* fall through */

			default:
				errno_exit("read");
			}
		}

		m_rgb=processFrame(buffers[0].start, buffers[0].length);
		break;

	case IO_METHOD_MMAP:
		CLEAR(buf);

		buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		buf.memory = V4L2_MEMORY_MMAP;

		if (-1 == xioctl(fd, VIDIOC_DQBUF, &buf)) {
			switch (errno) {
			case EAGAIN:
				return 0;

			case EIO:
				/* Could ignore EIO, see spec. */

				/* fall through */

			default:
				errno_exit("VIDIOC_DQBUF");
			}
		}

		assert(buf.index < n_buffers);

		m_rgb=processFrame(buffers[buf.index].start, buf.bytesused);

		if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
			errno_exit("VIDIOC_QBUF");
		break;

	case IO_METHOD_USERPTR:
		CLEAR(buf);

		buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		buf.memory = V4L2_MEMORY_USERPTR;

		if (-1 == xioctl(fd, VIDIOC_DQBUF, &buf)) {
			switch (errno) {
			case EAGAIN:
				return 0;

			case EIO:
				/* Could ignore EIO, see spec. */

				/* fall through */

			default:
				errno_exit("VIDIOC_DQBUF");
			}
		}

		for (i = 0; i < n_buffers; ++i)
			if (buf.m.userptr == (unsigned long)buffers[i].start
			    && buf.length == buffers[i].length)
				break;

		assert(i < n_buffers);

		m_rgb = processFrame((void *)buf.m.userptr, buf.bytesused);

		if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
			errno_exit("VIDIOC_QBUF");
		break;
	}

	return 1;
}

static void captureVideo(void)
{
	unsigned int count;
	char ch;
	pr_debug(DSP_INFO,"%s:\n", __func__);

	pr_debug(DSP_DEBUG,"win_size=%ds, fps=%d\n", getSampleWindow(),getFPS());

	initCQs();
	count = frame_count?frame_count:0xffffffff;
	while (count-- > 0) {
		for (;;) {
			fd_set fds;
			struct timeval tv;
			int r;

			FD_ZERO(&fds);
			FD_SET(fd, &fds);

			/* Timeout. */
			tv.tv_sec = 2;
			tv.tv_usec = 0;

			r = select(fd + 1, &fds, NULL, NULL, &tv);

			if (-1 == r) {
				if (EINTR == errno)
					continue;
				errno_exit("select");
			}

			if (0 == r) {
				fprintf(stderr, "select timeout\n");
				exit(EXIT_FAILURE);
			}

			if( (ch=cvWaitKey(1)) =='q') //this waitkey pause can make CV display visible
				goto exit;

			if (readFrame())
				break;
			/* EAGAIN - continue select loop. */
		}
	}
exit:;
	deinitCQs();
}

static void stopCapturing(void)
{
	enum v4l2_buf_type type;

	pr_debug(DSP_INFO,"%s:\n", __func__);

	switch (io) {
	case IO_METHOD_READ:
		/* Nothing to do. */
		break;

	case IO_METHOD_MMAP:
	case IO_METHOD_USERPTR:
		type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		if (-1 == xioctl(fd, VIDIOC_STREAMOFF, &type))
			errno_exit("VIDIOC_STREAMOFF");
		break;
	}
}

static void startCapturing(void)
{
	unsigned int i;
	enum v4l2_buf_type type;
	int err;

	pr_debug(DSP_INFO,"%s:\n", __func__);
	pr_debug(DSP_INFO,"\tn_buffers: %d\n", n_buffers);

	switch (io) {
	case IO_METHOD_READ:
		/* Nothing to do. */
		break;

	case IO_METHOD_MMAP:
		for (i = 0; i < n_buffers; ++i) {
			struct v4l2_buffer buf;

			pr_debug(DSP_INFO,"\ti: %d\n", i);

			CLEAR(buf);
			buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
			buf.memory = V4L2_MEMORY_MMAP;
			buf.index = i;

			pr_debug(DSP_INFO,"\tbuf.index: %d\n", buf.index);

			err == xioctl(fd, VIDIOC_QBUF, &buf);
			pr_debug(DSP_INFO,"\terr: %d\n", err);

			if (-1 == err)
				errno_exit("VIDIOC_QBUF");

			pr_debug(DSP_INFO,"\tbuffer queued!\n");
		}

		pr_debug(DSP_INFO,"Before STREAMON\n");
		type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		if (-1 == xioctl(fd, VIDIOC_STREAMON, &type))
			errno_exit("VIDIOC_STREAMON");
		pr_debug(DSP_INFO,"After STREAMON\n");
		break;

	case IO_METHOD_USERPTR:
		for (i = 0; i < n_buffers; ++i) {
			struct v4l2_buffer buf;

			CLEAR(buf);
			buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
			buf.memory = V4L2_MEMORY_USERPTR;
			buf.index = i;
			buf.m.userptr = (unsigned long)buffers[i].start;
			buf.length = buffers[i].length;

			if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
				errno_exit("VIDIOC_QBUF");
		}
		type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		if (-1 == xioctl(fd, VIDIOC_STREAMON, &type))
			errno_exit("VIDIOC_STREAMON");
		break;
	}
}

static void uninit_device(void)
{
	unsigned int i;

	pr_debug(DSP_INFO,"%s:\n", __func__);

	switch (io) {
	case IO_METHOD_READ:
		free(buffers[0].start);
		break;

	case IO_METHOD_MMAP:
		for (i = 0; i < n_buffers; ++i)
			if (-1 == munmap(buffers[i].start, buffers[i].length))
				errno_exit("munmap");
		break;

	case IO_METHOD_USERPTR:
		for (i = 0; i < n_buffers; ++i)
			free(buffers[i].start);
		break;
	}

	free(buffers);
}

static void init_read(unsigned int buffer_size)
{
	pr_debug(DSP_INFO,"%s:\n", __func__);

	buffers = calloc(1, sizeof(*buffers));

	if (!buffers) {
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}

	buffers[0].length = buffer_size;
	buffers[0].start = malloc(buffer_size);

	if (!buffers[0].start) {
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}
}

static void init_mmap(void)
{
	struct v4l2_requestbuffers req;

	pr_debug(DSP_INFO,"%s:\n", __func__);

	CLEAR(req);

	req.count = 4;
	req.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	req.memory = V4L2_MEMORY_MMAP;

	if (-1 == xioctl(fd, VIDIOC_REQBUFS, &req)) {
		if (EINVAL == errno) {
			fprintf(stderr, "%s does not support "
				 "memory mapping\n", dev_name);
			exit(EXIT_FAILURE);
		} else {
			errno_exit("VIDIOC_REQBUFS");
		}
	}
	pr_debug(DSP_INFO,"\treq.count: %d\n", req.count);
	pr_debug(DSP_INFO,"\treq.type: %d\n", req.type);
	pr_debug(DSP_INFO,"\treq.memory: %d\n", req.memory);
	pr_debug(DSP_INFO,"\n");

	if (req.count < 2) {
		fprintf(stderr, "Insufficient buffer memory on %s\n",
			 dev_name);
		exit(EXIT_FAILURE);
	}

	buffers = calloc(req.count, sizeof(*buffers));

	if (!buffers) {
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}

	for (n_buffers = 0; n_buffers < req.count; ++n_buffers) {
		struct v4l2_buffer buf;

		CLEAR(buf);

		buf.type        = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		buf.memory      = V4L2_MEMORY_MMAP;
		buf.index       = n_buffers;

		if (-1 == xioctl(fd, VIDIOC_QUERYBUF, &buf))
			errno_exit("VIDIOC_QUERYBUF");

		pr_debug(DSP_INFO,"\tbuf.index: %d\n", buf.index);
		pr_debug(DSP_INFO,"\tbuf.type: %d\n", buf.type);
		pr_debug(DSP_INFO,"\tbuf.bytesused: %d\n", buf.bytesused);
		pr_debug(DSP_INFO,"\tbuf.flags: %d\n", buf.flags);
		pr_debug(DSP_INFO,"\tbuf.field: %d\n", buf.field);
		pr_debug(DSP_INFO,"\tbuf.timestamp.tv_sec: %ld\n", (long) buf.timestamp.tv_sec);
		pr_debug(DSP_INFO,"\tbuf.timestamp.tv_usec: %ld\n", (long) buf.timestamp.tv_usec);
		pr_debug(DSP_INFO,"\tbuf.timecode.type: %d\n", buf.timecode.type);
		pr_debug(DSP_INFO,"\tbuf.timecode.flags: %d\n", buf.timecode.flags);
		pr_debug(DSP_INFO,"\tbuf.timecode.frames: %d\n", buf.timecode.frames);
		pr_debug(DSP_INFO,"\tbuf.timecode.seconds: %d\n", buf.timecode.seconds);
		pr_debug(DSP_INFO,"\tbuf.timecode.minutes: %d\n", buf.timecode.minutes);
		pr_debug(DSP_INFO,"\tbuf.timecode.hours: %d\n", buf.timecode.hours);
		pr_debug(DSP_INFO,"\tbuf.timecode.userbits: %d,%d,%d,%d\n",
				buf.timecode.userbits[0],
				buf.timecode.userbits[1],
				buf.timecode.userbits[2],
				buf.timecode.userbits[3]);
		pr_debug(DSP_INFO,"\tbuf.sequence: %d\n", buf.sequence);
		pr_debug(DSP_INFO,"\tbuf.memory: %d\n", buf.memory);
		pr_debug(DSP_INFO,"\tbuf.m.offset: %d\n", buf.m.offset);
		pr_debug(DSP_INFO,"\tbuf.length: %d\n", buf.length);
		pr_debug(DSP_INFO,"\tbuf.input: %d\n", buf.input);
		pr_debug(DSP_INFO,"\n");

		buffers[n_buffers].length = buf.length;
		buffers[n_buffers].start =
			mmap(NULL /* start anywhere */,
			      buf.length,
			      PROT_READ | PROT_WRITE /* required */,
			      MAP_SHARED /* recommended */,
			      fd, buf.m.offset);

		if (MAP_FAILED == buffers[n_buffers].start)
			errno_exit("mmap");
	}
}

static void init_userp(unsigned int buffer_size)
{
	struct v4l2_requestbuffers req;
	pr_debug(DSP_INFO,"%s:\n", __func__);

	CLEAR(req);

	req.count  = 4;
	req.type   = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	req.memory = V4L2_MEMORY_USERPTR;

	if (-1 == xioctl(fd, VIDIOC_REQBUFS, &req)) {
		if (EINVAL == errno) {
			fprintf(stderr, "%s does not support "
				 "user pointer i/o\n", dev_name);
			exit(EXIT_FAILURE);
		} else {
			errno_exit("VIDIOC_REQBUFS");
		}
	}

	buffers = calloc(4, sizeof(*buffers));

	if (!buffers) {
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}

	for (n_buffers = 0; n_buffers < 4; ++n_buffers) {
		buffers[n_buffers].length = buffer_size;
		buffers[n_buffers].start = malloc(buffer_size);

		if (!buffers[n_buffers].start) {
			fprintf(stderr, "Out of memory\n");
			exit(EXIT_FAILURE);
		}
	}
}

static void init_device(void)
{
	struct v4l2_capability cap;
	struct v4l2_cropcap cropcap;
	struct v4l2_crop crop;
	struct v4l2_format fmt;
	unsigned int min;

	pr_debug(DSP_INFO,"%s:\n", __func__);

	if (-1 == xioctl(fd, VIDIOC_QUERYCAP, &cap)) {
		if (EINVAL == errno) {
			fprintf(stderr, "%s is no V4L2 device\n",
				 dev_name);
			exit(EXIT_FAILURE);
		} else {
			errno_exit("VIDIOC_QUERYCAP");
		}
	}

	pr_debug(DSP_INFO,"\tdriver: %s\n"
		 "\tcard: %s \n"
		 "\tbus_info: %s\n",
			cap.driver, cap.card, cap.bus_info);
	pr_debug(DSP_INFO,"\tversion: %u.%u.%u\n",
			(cap.version >> 16) & 0xFF,
			(cap.version >> 8) & 0xFF,
			cap.version & 0xFF);
	pr_debug(DSP_INFO,"\tcapabilities: 0x%08x\n", cap.capabilities);

	if (!(cap.capabilities & V4L2_CAP_VIDEO_CAPTURE)) {
		fprintf(stderr, "%s is no video capture device\n",
			 dev_name);
		exit(EXIT_FAILURE);
	}

	switch (io) {
	case IO_METHOD_READ:
		if (!(cap.capabilities & V4L2_CAP_READWRITE)) {
			fprintf(stderr, "%s does not support read i/o\n",
				 dev_name);
			exit(EXIT_FAILURE);
		}
		break;

	case IO_METHOD_MMAP:
	case IO_METHOD_USERPTR:
		if (!(cap.capabilities & V4L2_CAP_STREAMING)) {
			fprintf(stderr, "%s does not support streaming i/o\n",
				 dev_name);
			exit(EXIT_FAILURE);
		}
		break;
	}


	/* Select video input, video standard and tune here. */


	CLEAR(cropcap);

	cropcap.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;

	if (0 == xioctl(fd, VIDIOC_CROPCAP, &cropcap)) {
		pr_debug(DSP_INFO,"\tcropcap.type: %d\n", cropcap.type);
		pr_debug(DSP_INFO,"\tcropcap.bounds.left: %d\n", cropcap.bounds.left);
		pr_debug(DSP_INFO,"\tcropcap.bounds.top: %d\n", cropcap.bounds.top);
		pr_debug(DSP_INFO,"\tcropcap.bounds.width: %d\n", cropcap.bounds.width);
		pr_debug(DSP_INFO,"\tcropcap.bounds.height: %d\n", cropcap.bounds.height);

		pr_debug(DSP_INFO,"\tcropcap.defrect.left: %d\n", cropcap.defrect.left);
		pr_debug(DSP_INFO,"\tcropcap.defrect.top: %d\n", cropcap.defrect.top);
		pr_debug(DSP_INFO,"\tcropcap.defrect.width: %d\n", cropcap.defrect.width);
		pr_debug(DSP_INFO,"\tcropcap.defrect.height: %d\n", cropcap.defrect.height);

		pr_debug(DSP_INFO,"\tcropcap.pixelaspect.numerator: %d\n", cropcap.pixelaspect.numerator);
		pr_debug(DSP_INFO,"\tcropcap.pixelaspect.denominator: %d\n", cropcap.pixelaspect.denominator);
		pr_debug(DSP_INFO,"\n");

		CLEAR(crop);
		crop.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		crop.c = cropcap.defrect; /* reset to default */

		if (-1 == xioctl(fd, VIDIOC_S_CROP, &crop)) {
			switch (errno) {
			case EINVAL:
				/* Cropping not supported. */
				break;
			default:
				/* Errors ignored. */
				pr_debug(DSP_INFO,"\tcropping not supported\n");
				break;
			}
		}
	} else {
		/* Errors ignored. */
	}


	CLEAR(fmt);

	fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	if (force_format) {
		fmt.fmt.pix.width       = FORCED_WIDTH;
		fmt.fmt.pix.height      = FORCED_HEIGHT;
		fmt.fmt.pix.pixelformat = FORCED_FORMAT;
		fmt.fmt.pix.field       = FORCED_FIELD;

		pr_debug(DSP_INFO,"\tfmt.fmt.pix.pixelformat: %c,%c,%c,%c\n",
				fmt.fmt.pix.pixelformat & 0xFF,
				(fmt.fmt.pix.pixelformat >> 8) & 0xFF,
				(fmt.fmt.pix.pixelformat >> 16) & 0xFF,
				(fmt.fmt.pix.pixelformat >> 24) & 0xFF
				);
		pr_debug(DSP_INFO,"\n");
		/*
		This, VIDIOC_S_FMT, is one time setting, so you can NOT set twice.
		You can't change the format while buffers are allocated.
		You need to free the buffers before, using VIDIOC_REQBUFS with a buffer
		count of 0.
		*/
		if (-1 == xioctl(fd, VIDIOC_S_FMT, &fmt))
			errno_exit("VIDIOC_S_FMT");

		/* Note VIDIOC_S_FMT may change width and height. */
	} else {
		/* Preserve original settings as set by v4l2-ctl for example */
		if (-1 == xioctl(fd, VIDIOC_G_FMT, &fmt))
			errno_exit("VIDIOC_G_FMT");

		fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		pr_debug(DSP_INFO,"\tfmt.fmt.pix.pixelformat: %c,%c,%c,%c\n",
				fmt.fmt.pix.pixelformat & 0xFF,
				(fmt.fmt.pix.pixelformat >> 8) & 0xFF,
				(fmt.fmt.pix.pixelformat >> 16) & 0xFF,
				(fmt.fmt.pix.pixelformat >> 24) & 0xFF
				);

	}

	/* Buggy driver paranoia. */
	min = fmt.fmt.pix.width * 2;
	if (fmt.fmt.pix.bytesperline < min)
		fmt.fmt.pix.bytesperline = min;
	min = fmt.fmt.pix.bytesperline * fmt.fmt.pix.height;
	if (fmt.fmt.pix.sizeimage < min)
		fmt.fmt.pix.sizeimage = min;

	extra_cam_setting(fd);

	switch (io) {
	case IO_METHOD_READ:
		init_read(fmt.fmt.pix.sizeimage);
		break;

	case IO_METHOD_MMAP:
		init_mmap();
		break;

	case IO_METHOD_USERPTR:
		init_userp(fmt.fmt.pix.sizeimage);
		break;
	}
}

static void close_device(void)
{
	pr_debug(DSP_INFO,"%s:\n", __func__);

	if (-1 == close(fd))
		errno_exit("close");

	fd = -1;
}

static void open_device(void)
{
	struct stat st;

	pr_debug(DSP_INFO,"+%s:\n", __func__);

	if (-1 == stat(dev_name, &st)) {
		fprintf(stderr, "Cannot identify '%s': %d, %s\n",
			 dev_name, errno, strerror(errno));
		exit(EXIT_FAILURE);
	}

	if (!S_ISCHR(st.st_mode)) {
		fprintf(stderr, "%s is no device\n", dev_name);
		exit(EXIT_FAILURE);
	}

	fd = open(dev_name, O_RDWR /* required */ | O_NONBLOCK, 0);

	if (-1 == fd) {
		fprintf(stderr, "Cannot open '%s': %d, %s\n",
			 dev_name, errno, strerror(errno));
		exit(EXIT_FAILURE);
	}
}

static void usage(FILE *fp, int argc, char **argv)
{
	fprintf(fp,
		 "Usage: %s [options]\n\n"
		 "Version 0.6\n"
		 "Options:\n"
		 "-h | --help          Print this message\n"
		 "-c | --count         Number of frames to grab [%i]\n"
		 "-d | --device name   Video device name [%s]\n"
		 "-e | --errorlevel    error level [%i]\n"
		 "-i | --samplefile    raw trace sample file name\n"
		 "-p | --mmap          Use memory mapped buffers [default]\n"
		 "-m | --min           minimum sample period in seconds[%i]\n"
		 "-M | --max           max sample period in seconds[%i]\n"
		 "-r | --read          Use read() calls\n"
		 "-R | --RADICAL       ICA-RADICAL on/off[%d]\n"
		 "-u | --userp         Use application allocated buffers\n"
		 "-O | --stdout        Outputs stream to stdout\n"
		 "-o | --output        debug file stream to stdout\n"
		 "-w | --window        working window in second[%i]\n"
		 "-f | --format        Force format to 640x480 YUYV\n"
		 "-F | --fps           fps [%i]\n"
		 "-s | --steps         steps in seconds[%i]\n"
		 "-v | --verbose       Verbose output\n"
		 "",
		 argv[0], frame_count, dev_name, getErrorLevel(), getMinSampleTime(),
		getMaxSampleTime(), getRADICAL(), getSampleWindow(),getFPS(),getStepping());
}

static const char short_options[] = "c:d:e:f:F:hi:m:M:o:prR:s:uvw:";

static const struct option
long_options[] = {
	{ "help",   no_argument,       NULL, 'h' },
	{ "count",  required_argument, NULL, 'c' },
	{ "device", required_argument, NULL, 'd' },
	{ "errorlevel", required_argument, NULL, 'e' },
	{ "format", no_argument,       NULL, 'f' },
	{ "fps",  		required_argument, NULL, 'F' },
	{ "sample-file",required_argument, NULL, 'i' },
	{ "min",		required_argument, NULL, 'm' },
	{ "max",		required_argument, NULL, 'M' },
	{ "output", 	required_argument, NULL, 'o' },
	{ "stdout", 	required_argument, NULL, 'O' },
	{ "mmap",   no_argument,       NULL, 'p' },
	{ "read",   no_argument,       NULL, 'r' },
	{ "RADICAL", 	required_argument, NULL, 'R' },
	{ "step",  		required_argument, NULL, 's' },
	{ "userp",  no_argument,       NULL, 'u' },
	{ "verbose", 	no_argument,       NULL, 'v' },
	{ "window",		required_argument, NULL, 'w' },
	{ 0, 0, 0, 0 }
};

static int option(int argc, char **argv)
{
	int r=0;
	printf("+%s:\n",__func__);
	dev_name = "/dev/video0";

	for (;;) {
		int idx;
		int c;

		c = getopt_long(argc, argv,
				short_options, long_options, &idx);

		if (-1 == c)
			break;

		switch (c) {
		case 0: /* getopt_long() flag */
			break;
		case 'c':
			errno = 0;
			frame_count = strtol(optarg, NULL, 0);
			if (errno)
				errno_exit(optarg);
			break;

		case 'd':
			dev_name = optarg;
			break;

		case 'e':
			errno = 0;
			setErrorLevel(strtol(optarg, NULL, 0));
			if (errno)
				errno_exit(optarg);
			break;

		case 'f':
			force_format++;
			break;

		case 'F':
			errno = 0;
			setFPS(strtol(optarg, NULL, 0));
			if (errno)
				errno_exit(optarg);
			break;

		case 'h':
			usage(stdout, argc, argv);
			r=-1;
			break;

		case 'o':
			out_buf++;
			break;

		case 'p':
			io = IO_METHOD_MMAP;
			break;

		case 'r':
			io = IO_METHOD_READ;
			break;

		case 'R':
			errno = 0;
			setRADICAL(strtol(optarg, NULL, 0));
			printf("DSP RAIDCAL [%d].\n",getRADICAL());
			if (errno)
				errno_exit(optarg);
			break;

		case 'u':
			io = IO_METHOD_USERPTR;
			break;

		case 'v':
			setVerbose(1);
			break;

		case 'w':
			errno = 0;
			setSampleWindow(strtol(optarg, NULL, 0));
			if (errno)
				errno_exit(optarg);
			break;

		default:
			usage(stderr, argc, argv);
			r=-1;
		}
	}
	return r;
}

int main(int argc, char **argv)
{
	if(option(argc,argv)<0){
		pr_debug(DSP_DEBUG, "Wrong args!!!\n");
		exit(EXIT_FAILURE);
	}
	if(_init_dsp()){
		pr_debug(DSP_DEBUG, "DSP init failed!!!\n");
		exit(EXIT_FAILURE);
	}
	//start_logging("mylog.txt");
	open_device();
	init_device();

	cvNamedWindow(windowname,CV_WINDOW_AUTOSIZE);

	startCapturing();
	captureVideo();
	stopCapturing();

	cvDestroyWindow(windowname);

	uninit_device();
	close_device();
	_deinit_dsp();

	//fprintf(stderr, "\n");
	return 0;
}
