#include <stdio.h>

#include <cv.h>
#include <highgui.h>

extern void process_webcam_frames();
extern int minsize;
extern int maxsize;

int main(int argc, char* argv[])
{
	IplImage* img = 0;

	if(argc==1)
	{
		printf("Copyright (c) 2013, Nenad Markus\n");
		printf("All rights reserved.\n\n");

		process_webcam_frames();
	}
	else if(argc==2)
	{
		printf("Copyright (c) 2013, Nenad Markus\n");
		printf("All rights reserved.\n\n");

		sscanf(argv[1], "%d", &minsize);

		process_webcam_frames();
	}
	else if(argc==3)
	{
		sscanf(argv[1], "%d", &minsize);

		img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);
		if(!img)
		{
			printf("* cannot load image!\n");
			return 1;
		}

		process_image(img, 0, 1);

		cvReleaseImage(&img);
	}
	else if(argc==4)
	{
		sscanf(argv[1], "%d", &minsize);

		img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);
		if(!img)
		{
			printf("* cannot load image!\n");
			return 1;
		}

		process_image(img, 1, 0);

		//
		cvSaveImage(argv[3], img, 0);

		//
		cvReleaseImage(&img);
	}
	else
		return 1;

	return 0;
}
