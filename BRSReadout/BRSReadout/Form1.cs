using System;
using System.Collections;
using System.Collections.Generic;
using System.Configuration;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Windows.Forms;
using TwinCAT.Ads;


namespace BRSReadout
{
    /*
     * 
     * The primary declaration of the Form1 class that handels the primary functions and structure of the program
     * 
     */
    public partial class Form1 : Form
    {
        public bool twinCatBool = ("true" == ConfigurationManager.AppSettings.Get("twinCat"));
        public string cameraType = ConfigurationManager.AppSettings.Get("camera");

        public Form2 graphWindow = new Form2();

        public const int datalen = 4096;
        public const int datamax = 4096;
        public bool firstFrame = true;
        public int firstValueCounter = 0;
        public double zeroValue = 0;
        public double refZeroValue = 0;
        public double refValue = 0;

        public int frameCount = 0;

        Camera myCamera;
        DataConsumerDelegate[] consumerd;

        Stopwatch myStopwatch;
        DataWriting myDataWriter;

        public string toWrite;
        private RawData currentData;

        int startIndex2 = 1800; //Beginning of right pattern
        int startIndex2Ref = 0; //Beginning of right pattern
        bool noMoreData = false;
        volatile bool gquitting = false;
        public volatile bool gRecord = false;
        public volatile bool gDamp = false;

        double gsampfreq = Math.Pow(10,6)/double.Parse(ConfigurationManager.AppSettings.Get("cameraExposureTime"))/2;
        double gmastersampfreq = Math.Pow(10, 6) / double.Parse(ConfigurationManager.AppSettings.Get("cameraExposureTime"))/2;

        public int gMaiIx = 1700;
        public volatile int Frameco = 0;
        public volatile int gNoPeaks = 1, gOldNoPeaks = 0;
        int Trig = 0;
        int gsampfrequcounter = 0;
        int gsampfrequcountermax = 1000;
        int gsampfrequcountermaxlater = 1000;

        ushort[] refFrame = new ushort[4096];


        public static ushort[] frame = new ushort[4096];
        ushort[] refFrameRef = new ushort[4096];
        volatile bool dataWritingThreadOn;
        Thread dataWritingThread;
        Thread cameraThread;
        private Queue<PeakQueueItem> dataWritingQueue;
        public volatile SyncEvents dataWritingSyncEvent;

        public static object graphLock = new object();
        public static Queue<graphData> graphQueue = new Queue<graphData>();
        public static AutoResetEvent graphSignal = new AutoResetEvent(false);
        Thread graphThread;
        public double graphSum;

        DateTime DTFrameCo0;
        double dayFrameCo0;

        double voltagewrite = 0;

        static int gFrames = 10; //Amount of frames collected before fitting and downsampling 
        public static double[,] newdata = new double[gFrames, 2];
        static TcAdsClient tcAds = new TcAdsClient();
        static AdsStream ds = new AdsStream(16);

        double x = 0;
        double xSquar = 0;
        double xCube = 0;
        double xFourth = 0;

        double[] yHighPass = new double[3];
        double[] xHighPass = new double[3];
        double[] yBandHighPass = new double[3];
        double[] xBandHighPass = new double[3];
        double[] yBandLowPass = new double[3];
        double[] xBandLowPass = new double[3];
        double[] highCoeff = new double[5];
        double[] bandLowCoeff = new double[5];
        double[] bandHighCoeff = new double[5];

        double[] yLowPass1 = new double[3];
        double[] xLowPass1 = new double[3];
        double[] yLowPass2 = new double[3];
        double[] xLowPass2 = new double[3];
        double[] lowCoeff = new double[5];


        int lightSourceStatus;
        int cameraStatus = 1;

        int displayCount = 0;

        public struct graphData
        {
            public ushort[] frame;
            public double angle;
            public graphData(ushort[] inFrame, double inAngle)
            {
                angle = inAngle;
                frame = inFrame;
            }
        }

        public Form1()
        {
            try
            {
                //Calls calculation method for filters
                highCoeff = filterCoeff(0.0005, 2000.0 / gFrames, "High");
                lowCoeff = filterCoeff(0.1, 2000.0 / gFrames, "Low");
                bandHighCoeff = filterCoeff(0.002, 2000.0 / gFrames, "High");
                bandLowCoeff = filterCoeff(0.02, 2000.0 / gFrames, "Low");//2*10^-2

                InitializeComponent(); // Initializes the visual components
                SetSize(); // Sets up the size of the window and the corresponding location of the components
                Frameco = 0;
                DTFrameCo0 = DateTime.Now;
                dayFrameCo0 = DayNr.GetDayNr(DTFrameCo0);

                toWrite = "";
                currentData = new RawData(gFrames);
                //Initialization of TwinCAT connection. 851 is the port for PLC runtime 1
                //tcAds.Dispose();
                if (twinCatBool)
                {
                    tcAds.Connect(851);
                }

                myStopwatch = new Stopwatch();
                dataWritingSyncEvent = new SyncEvents();

                consumerd = new DataConsumerDelegate[2];
                consumerd[0] = new DataConsumerDelegate(new DataConsumerDelegate.outputDelegate(showStatistics), true);
                consumerd[1] = new DataConsumerDelegate(new DataConsumerDelegate.outputDelegate(Pattern), true);
                consumerd[1].myThread.Priority = ThreadPriority.Highest;
                for (int i = 0; i < consumerd.Length; i++)
                {
                    consumerd[i].myThread.Start();
                }
                dataWritingQueue = new Queue<PeakQueueItem>();

                dataWritingThreadOn = true;
                dataWritingThread = new Thread(dataWrite);
                dataWritingThread.Priority = ThreadPriority.Highest;
                dataWritingThread.Start();
                Thread.Sleep(10);

                cameraThread = new Thread(initCamera);
                cameraThread.Priority = ThreadPriority.AboveNormal;
                cameraThread.Start();
                string curDirec = System.IO.Path.GetDirectoryName(Application.ExecutablePath);
                curDirec = curDirec.Replace("\\bin\\Debug", "");
                System.Diagnostics.Process.Start(curDirec + "\\AffinitySet.bat");
            }
            catch (Exception ex)
            {
                EmailError.emailAlert(ex);
                throw (ex);
            }
        }
        // Local data writing
        void dataWrite()
        {
            int f;
            double diff;
            double reff;
            StringBuilder ds1, ds2;
            PeakQueueItem curItem;
            double[,] dataA;
            long[] timeStampA;
            long curTimeStamp;


            double[] refLP;
            double[,] values;

            values = new double[10000, 2];
            curItem = new PeakQueueItem();
            while (dataWritingThreadOn)
            {
                if (gquitting)
                {
                    return;
                }

                while ((WaitHandle.WaitAny(dataWritingSyncEvent.EventArray) != 1))
                {
                    if (gquitting)
                    {
                        return;
                    }

                    lock (((ICollection)dataWritingQueue).SyncRoot)
                    {
                        curItem = dataWritingQueue.Dequeue();
                        if (dataWritingQueue.Count >= 1)
                        {
                            dataWritingSyncEvent.NewItemEvent.Set();
                        }

                    }
                    ds1 = new StringBuilder(100000);
                    ds2 = new StringBuilder(100000);
                    timeStampA = curItem.TimeStamp;
                    dataA = curItem.Data;
                    curTimeStamp = 0;
                    for (f = 0; f < gFrames; f++)
                    {
                        curTimeStamp = timeStampA[f];
                        ds1.AppendFormat("{0,10:F2}", curTimeStamp.ToString());
                        ds2.AppendFormat("{0,10:F2}", curTimeStamp.ToString());
                        diff = dataA[f, 0];
                        reff = dataA[f, 1];
                        values[f, 0] = diff;
                        values[f, 1] = reff;
                    }
                    refLP = lp(values, gFrames);
                    dataSend(refLP);

                    f = f - 1;
                    if (gRecord)
                    {
                        myDataWriter.Write(curTimeStamp / gmastersampfreq / 3600 / 24 + dayFrameCo0, refLP);
                    }
                    if (Application.OpenForms.OfType<Form2>().Count() == 1)
                    {
                        for (int i = 0; i < newdata.GetLength(0); i++)
                        {
                            if (newdata[i, 0] > 0)
                            {
                                graphSum += newdata[i, 0];
                            }
                        }
                        graphSum = graphSum / ((double)newdata.GetLength(0));
                        graphData outData = new graphData(frame, graphSum);
                        lock (graphLock)
                        {
                            graphQueue.Enqueue(outData);
                        }
                        graphSignal.Set();
                        frameCount = 0;
                        graphSum = 0;
                    }
                    else
                    {
                        frameCount++;
                    }


                }

            }
        }
        /*
         *  Takes a two dimensional array of doubles and int and returns an array of doubles
         *  Low pass filter for data
         */
        double[] lp(double[,] data, int frames)
        {
            double[] output = new double[gFrames];
            double mainSum = 0;
            double refSum = 0;
            for (int i = 0; i < gFrames; i++)
            {
                xLowPass1[0] = data[i, 0];

                yLowPass1[0] = lowCoeff[3] * yLowPass1[1] + lowCoeff[4] * yLowPass1[2] + lowCoeff[0] * xLowPass1[0] + lowCoeff[1] * xLowPass1[1] + lowCoeff[2] * xLowPass1[2];

                xLowPass1[2] = xLowPass1[1];
                xLowPass1[1] = xLowPass1[0];
                yLowPass1[2] = yLowPass1[1];
                yLowPass1[1] = yLowPass1[0];

                xLowPass2[0] = data[i, 1];

                yLowPass2[0] = lowCoeff[3] * yLowPass2[1] + lowCoeff[4] * yLowPass2[2] + lowCoeff[0] * xLowPass2[0] + lowCoeff[1] * xLowPass2[1] + lowCoeff[2] * xLowPass2[2];

                xLowPass2[2] = xLowPass2[1];
                xLowPass2[1] = xLowPass2[0];
                yLowPass2[2] = yLowPass2[1];
                yLowPass2[1] = yLowPass2[0];

                mainSum += yLowPass1[0];
                refSum += yLowPass2[0];
            }
            output[0] = mainSum / gFrames;
            output[1] = refSum / gFrames;
            return output;
        }
        // This method takes new data from the camera and inserts it, with proper timing, into the queue to be processed by the consumerd
        // the data consumer delegate that processes and displays the data
        public void onNewData(ushort[] data)
        {
            int i;
            int curFrame;
            Frameco++;
            if (gquitting)
            {
                return;
            }

            gsampfrequcounter++;
            if (gsampfrequcounter == 1)
            {
                myStopwatch.Reset();
                myStopwatch.Start();
            }
            if (gsampfrequcounter == gsampfrequcountermax)
            {
                myStopwatch.Stop();
                TimeSpan ts = myStopwatch.Elapsed;
                gsampfreq = ((gsampfrequcountermax - 1.0) / ts.TotalMilliseconds * 1000.0);
                gsampfrequcounter = 0;
                if (gsampfrequcountermax != gsampfrequcountermaxlater)
                {
                    gsampfrequcountermax = gsampfrequcountermaxlater;
                }
            }
            // Enque the Data
            if (noMoreData)
            {
                return;
            }
            curFrame = currentData.AddData(data, Frameco); // Returns the number of frames in the RawData object

            if (curFrame == gFrames)
            {
                try
                {
                    for (i = 0; i < consumerd.Length; i++)
                    {
                        lock (((ICollection)consumerd[i].myQueue).SyncRoot)
                        {
                            consumerd[i].myQueue.Enqueue(currentData);
                        }
                        consumerd[i].mySyncEvent.NewItemEvent.Set();
                    }
                }
                catch (Exception ex)
                {
                    EmailError.emailAlert(ex);
                    throw (ex);
                }
                currentData = new RawData(gFrames);
            }
        }
        //Camera initialization
        public void initCamera()
        {
            try
            {
                myCamera = new Camera();
                myCamera.cameraInit(cameraType);
                Camera.dataDelegate dd = new Camera.dataDelegate(onNewData);
                myCamera.startFrameGrab(0x8888, 0, dd, cameraType);
                cameraStatus = 1;
            }
            catch (Exception ex)
            {
                cameraStatus = 0;
                //Camera status bit
                ds = new AdsStream(4);
                BinaryWriter bw = new BinaryWriter(ds);
                bw.Write(cameraStatus);
                if (twinCatBool)
                {
                    tcAds.Write(0x4020, 40, ds);
                }
                EmailError.emailAlert(ex);
                throw (ex);
            }
        }

        //GUI display of queue length, time, and capacitor voltage
        private void showStatistics(RawData data)
        {
            if (displayCount == 10)
            {
                double ti;
                ti = data.TimeStamp(0) * 1.0 / gmastersampfreq;
                setTextBox1(ti.ToString("F1"));
                setTextBox4(String.Format("{0:0.00}", (voltagewrite)));
                displayCount = 0;
            }
            displayCount++;

        }
        //Calculates filter coefficients for a second order Butterworth of passed type
        double[] filterCoeff(double cutFreq, double sampFreq, string type)
        {
            double a1, a2, a3, b1, b2, c;
            double[] output = new double[5];
            if (type.Equals("Low") == true)
            {
                c = Math.Pow(Math.E, -2 * Math.PI * cutFreq / (sampFreq));
                a1 = Math.Pow((1.0 - c), 2);
                a2 = 0;
                a3 = 0;
                b1 = 2 * c;
                b2 = -Math.Pow(c, 2);
            }
            else if (type.Equals("High") == true)
            {
                c = Math.Pow(Math.E, -2 * Math.PI * cutFreq / (sampFreq));
                a1 = Math.Pow((1.0 + c) / 2, 2);
                a2 = -2 * a1;
                a3 = a1;
                b1 = 2 * c;
                b2 = -Math.Pow(c, 2);
                ;
            }
            else
            {
                a1 = 1; a2 = 0; a3 = 0; b1 = 0; b2 = 0;
            }
            output[0] = a1; output[1] = a2; output[2] = a3; output[3] = b1; output[4] = b2;
            return output;
        }
        //Prepares data to be sent to TwinCAT software
        private void dataSend(double[] data)
        {
            double tilt;
            double drift = 0;
            double velocity = 0;
            double angle = 0;
            double refAng = 0;
            angle = data[0];
            refAng = 58.33 * 60 * (data[1] - refZeroValue);
            //Has a DC subtraction to help filters. Should be about the center of the signal.
            if (firstValueCounter < 40)
            {
                firstValueCounter++;
                zeroValue = data[0];
                refZeroValue = data[1];
            }
            xHighPass[0] = angle - zeroValue;
            xBandLowPass[0] = angle - zeroValue;

            //Drift signal calculation, just scaled signal
            drift = 60 * (angle - 2100 + 80);
            //Drift rail logic
            if (Math.Abs(drift) >= 32760)
            {
                drift = Math.Sign(drift) * 32760;
            }

            //Tilt signal calculations, high pass at 10^-3 Hz then scaling.
            yHighPass[0] = highCoeff[3] * yHighPass[1] + highCoeff[4] * yHighPass[2] + highCoeff[0] * xHighPass[0] + highCoeff[1] * xHighPass[1] + highCoeff[2] * xHighPass[2];
            tilt = 0.729 * 5000 * yHighPass[0];

            //Capacitor signal calculations, low passed at Hz then high passed at Hz then differentiated and scaled
            yBandLowPass[0] = bandLowCoeff[3] * yBandLowPass[1] + bandLowCoeff[4] * yBandLowPass[2] + bandLowCoeff[0] * xBandLowPass[0] + bandLowCoeff[1] * xBandLowPass[1] + bandLowCoeff[2] * xBandLowPass[2];
            xBandHighPass[0] = yBandLowPass[0];
            yBandHighPass[0] = bandHighCoeff[3] * yBandHighPass[1] + bandHighCoeff[4] * yBandHighPass[2] + bandHighCoeff[0] * xBandHighPass[0] + bandHighCoeff[1] * xBandHighPass[1] + bandHighCoeff[2] * xBandHighPass[2];
            velocity = -2000000 * Math.Round((double)200 / gFrames) * (yBandHighPass[0] - yBandHighPass[1]);

            //Contact potential check. If above contact potential, the force is approximately proportional to V^2. If below, to -V*Vcontact.
            if (Math.Abs(velocity) > 00)
            {
                if (velocity > 0)
                {
                    velocity = 50 * Math.Sqrt(velocity);
                }
                else
                {
                    velocity = -50 * Math.Sqrt(-velocity);
                }
            }
            else
            {
                velocity *= -1.5;
            }

            if (Math.Abs(velocity) >= 30760)
            {
                velocity = Math.Sign(velocity) * 30760;
            }
            if (twinCatBool)
            {
                ds = new AdsStream(28);
                BinaryWriter bw = new BinaryWriter(ds);
                //Tilt signal. TwinCAT variable tilt at MW0.
                bw.Write((int)tilt);
                //Drift signal. TwinCAT variable drift at MW1
                bw.Write((int)drift);
                //Velocity signal. TwinCAT variable cap at MW2
                bw.Write((int)velocity);
                voltagewrite = velocity / 3276;
                //Reference signal. TwinCAT variable ref at MW3
                bw.Write((int)refAng);
                //C# pulse. Sets TwinCAT variable cPulse=1 at MW4
                bw.Write((int)1);
                //Light source status bit at MW5
                bw.Write((int)lightSourceStatus);
                //Camera status bit at MW6
                bw.Write((int)cameraStatus);
                tcAds.Write(0x4020, 0, ds);
            }

            xHighPass[2] = xHighPass[1];
            xHighPass[1] = xHighPass[0];
            yHighPass[2] = yHighPass[1];
            yHighPass[1] = yHighPass[0];
            xBandLowPass[2] = xBandLowPass[1];
            xBandLowPass[1] = xBandLowPass[0];
            yBandLowPass[2] = yBandLowPass[1];
            yBandLowPass[1] = yBandLowPass[0];
            xBandHighPass[2] = xBandHighPass[1];
            xBandHighPass[1] = xBandHighPass[0];
            yBandHighPass[2] = yBandHighPass[1];
            yBandHighPass[1] = yBandHighPass[0];
        }
        // This functions processes a pattern, obtains the data and send it out to be written
        private void Pattern(RawData data)
        {
            int nf;
            PeakQueueItem quI;
            int ql;
            long[] timestamps;
            double fitLength = 15; //Amount of pixels to be used in fit of correlation  - changed on 11/08/2017 by Krishna
            int halflength = (int)Math.Floor(fitLength / 2) + 1;  // increased by 1 pixel to allow two fits
            int length = 1500; //Length of patterns
            double[] crossCor = new double[(int)fitLength + 2];   // increased by 2 pixels to allow two fits
            int startIndex1 = 60; //Beginning of left pattern
            int startIndex1Ref = 60; //Beginning of left pattern
            int pixshift = 1;  // Direction of shift required to estimate slope of fit correction


            float sum = 0;
            double[] offset = new double[1];
            double[] fit = new double[4];
            double mu = 0;
            double mu1 = 0;
            double mu2 = 0;
            double N = fitLength;
            double y = 0;
            double xy = 0;
            double xxy = 0;
            double b, c, D, Db, Dc;// a,b, and c are the values we solve for then derive mu,sigma, and A from them.

            if (gquitting)
            {
                return;
            }

            Trig++;
            nf = gFrames;
            timestamps = new long[nf];
            newdata = new double[nf, 2];
            try
            {
                for (int frameNo = 0; frameNo < gFrames; frameNo++)
                {
                    frame = data.getData(frameNo);
                    if (firstFrame == true)
                    {
                        refFrame = frame;
                        firstFrame = false;
                        for (int j = 1600; j < frame.Length; j++)
                        {
                            if (frame[j] > 1200)//&& j < 1600)
                            {
                                startIndex2 = j - 30;
                                break;
                            }
                        }
                        for (int j = 0; j < frame.Length; j++)
                        {
                            if (frame[j] > 1200)//&& j < 1600)
                            {
                                startIndex2Ref = j - 30;
                                break;
                            }
                        }
                        for (int j = 0; j < fitLength; j++)
                        {
                            x += j;
                            xSquar += j * j;
                            xCube += j * j * j;
                            xFourth += j * j * j * j;
                        }
                    }
                    //Finds beginning of pattern using a threshold
                    if (frameNo == 0)
                    {
                        for (int j = 1600; j < frame.Length; j++)
                        {
                            if (frame[j] > 1200)//&& j < 1600)
                            {
                                startIndex1 = j - 30;
                                break;
                            }

                        }
                        if (startIndex1 < 0)
                        {
                            startIndex1 = 0;
                        }
                    }
                    //Cuts length of pattern down if the pattern extends beyond the frame
                    while (length + startIndex1 + halflength + 1 >= frame.Length)
                    {
                        length = (int)Math.Round(length / 1.1);
                    }
                    //Calcualtes the crosscorrelation between the two patterns at shifts ; first time
                    for (int k = -halflength; k <= halflength; k++)
                    {
                        sum = 0;
                        for (int m = 0; m < length; m++)
                        {
                            if ((m + startIndex1 + k) > 0 && (m + startIndex2) > 0)
                            {
                                sum += frame[m + startIndex1 + k] * refFrame[m + startIndex2];
                            }
                        }
                        if (sum == 0)
                        {
                            cameraStatus = 0;
                        }
                        else
                        {
                            cameraStatus = 1;
                        }
                        crossCor[k + halflength] = sum;
                    }
                    //Sums x,x^2,x^3,x^4,ln(y),x ln(y),x^2 ln(y)
                    y = 0;
                    xy = 0;
                    xxy = 0;
                    for (int j = 0; j < fitLength; j++)
                    {
                        y += Math.Log(crossCor[j + 1]);
                        xy += j * Math.Log(crossCor[j + 1]);
                        xxy += j * j * Math.Log(crossCor[j + 1]);
                    }
                    //Solves system of equations using Cramer's rule
                    D = N * (xSquar * xFourth - xCube * xCube) - x * (x * xFourth - xCube * xSquar) + xSquar * (x * xCube - xSquar * xSquar);
                    //Da = y * (xSquar * xFourth - xCube * xCube) - x * (xy * xFourth - xCube * xxy) + xSquar * (xy * xCube - xSquar * xxy);
                    Db = N * (xy * xFourth - xCube * xxy) - y * (x * xFourth - xCube * xSquar) + xSquar * (x * xxy - xy * xSquar);
                    Dc = N * (xSquar * xxy - xy * xCube) - x * (x * xxy - xy * xSquar) + y * (x * xCube - xSquar * xSquar);
                    //a = Da / D;
                    b = Db / D;
                    c = Dc / D;

                    mu1 = -b / (2 * c);

                    // If fit-center is to left of center of crosscor pattern, shift cross-cor pattern by 1 pixel to right or vice versa
                    if (mu1 < (halflength - 1))
                    {
                        pixshift = -1;
                    }
                    else
                    {
                        pixshift = 1;
                    }

                    //Redo the fit
                    y = 0;
                    xy = 0;
                    xxy = 0;
                    for (int j = 0; j < fitLength; j++)
                    {
                        y += Math.Log(crossCor[j + 1 + pixshift]);
                        xy += j * Math.Log(crossCor[j + 1 + pixshift]);
                        xxy += j * j * Math.Log(crossCor[j + 1 + pixshift]);
                    }
                    D = N * (xSquar * xFourth - xCube * xCube) - x * (x * xFourth - xCube * xSquar) + xSquar * (x * xCube - xSquar * xSquar);
                    Db = N * (xy * xFourth - xCube * xxy) - y * (x * xFourth - xCube * xSquar) + xSquar * (x * xxy - xy * xSquar);
                    Dc = N * (xSquar * xxy - xy * xCube) - x * (x * xxy - xy * xSquar) + y * (x * xCube - xSquar * xSquar);
                    b = Db / D;
                    c = Dc / D;

                    mu2 = -b / (2 * c);


                    mu = (halflength - 1) - (mu1 - (halflength - 1)) * pixshift / (mu2 - mu1);

                    newdata[frameNo, 0] = mu + startIndex1;
                    timestamps[frameNo] = data.TimeStamp(frameNo);





                    ////////////////////////////Reference Pattern/////////////////////
                    y = 0;
                    xy = 0;
                    xxy = 0;
                    //Finds beginning of pattern using a threshold

                    // 11/08/17 - Changed code to online subtraction of Ref signal - Krishna
                    //if (frameCount >= 10)
                    //{
                    //    frameCount = 0;

                    for (int j = 0; j < frame.Length; j++)
                    {
                        if (frame[j] > 1200)//&& j < 1600)
                        {
                            startIndex1Ref = j - 30;
                            lightSourceStatus = 1;
                            break;
                        }

                        if (j == frame.Length - 1)
                        {
                            lightSourceStatus = 0;
                        }
                    }
                    if (startIndex1Ref < halflength)
                    {
                        startIndex1Ref = 0;
                    }

                    //Calcualtes the crosscorrelation between the two patterns at shifts 
                    for (int k = -halflength; k <= halflength; k++)
                    {
                        sum = 0;
                        for (int m = 0; m < length; m++)
                        {
                            if ((m + startIndex1Ref + k) > 0 && (m + startIndex2Ref) > 0)
                            {
                                sum += frame[m + startIndex1Ref + k] * refFrame[m + startIndex2Ref];
                            }
                        }
                        crossCor[k + halflength] = sum;
                    }
                    //Sums x,x^2,x^3,x^4,ln(y),x ln(y),x^2 ln(y)

                    for (int j = 0; j < fitLength; j++)
                    {
                        y += Math.Log(crossCor[j + 1]);
                        xy += j * Math.Log(crossCor[j + 1]);
                        xxy += j * j * Math.Log(crossCor[j + 1]);
                    }
                    //Solves system of equations using Cramer's rule
                    D = N * (xSquar * xFourth - xCube * xCube) - x * (x * xFourth - xCube * xSquar) + xSquar * (x * xCube - xSquar * xSquar);
                    //Da = y * (xSquar * xFourth - xCube * xCube) - x * (xy * xFourth - xCube * xxy) + xSquar * (xy * xCube - xSquar * xxy);
                    Db = N * (xy * xFourth - xCube * xxy) - y * (x * xFourth - xCube * xSquar) + xSquar * (x * xxy - xy * xSquar);
                    Dc = N * (xSquar * xxy - xy * xCube) - x * (x * xxy - xy * xSquar) + y * (x * xCube - xSquar * xSquar);
                    //a = Da / D;
                    b = Db / D;
                    c = Dc / D;

                    mu1 = -b / (2 * c);

                    // If fit-center is to left of center of crosscor pattern, shift cross-cor pattern by 1 pixel to right or vice versa
                    if (mu1 < (halflength - 1))
                    {
                        pixshift = -1;
                    }
                    else
                    {
                        pixshift = 1;
                    }

                    //Redo the fit
                    y = 0;
                    xy = 0;
                    xxy = 0;
                    for (int j = 0; j < fitLength; j++)
                    {
                        y += Math.Log(crossCor[j + 1 + pixshift]);
                        xy += j * Math.Log(crossCor[j + 1 + pixshift]);
                        xxy += j * j * Math.Log(crossCor[j + 1 + pixshift]);
                    }
                    D = N * (xSquar * xFourth - xCube * xCube) - x * (x * xFourth - xCube * xSquar) + xSquar * (x * xCube - xSquar * xSquar);
                    Db = N * (xy * xFourth - xCube * xxy) - y * (x * xFourth - xCube * xSquar) + xSquar * (x * xxy - xy * xSquar);
                    Dc = N * (xSquar * xxy - xy * xCube) - x * (x * xxy - xy * xSquar) + y * (x * xCube - xSquar * xSquar);
                    b = Db / D;
                    c = Dc / D;

                    mu2 = -b / (2 * c);

                    mu = (halflength - 1) - (mu1 - (halflength - 1)) * pixshift / (mu2 - mu1);

                    newdata[frameNo, 1] = mu + startIndex1Ref;
                    refValue = mu + startIndex1Ref;
                    newdata[frameNo, 0] = newdata[frameNo, 0] - refValue;

                }
            }
            catch (Exception ex)
            {
                EmailError.emailAlert(ex);
                throw (ex);
            }

            quI = new PeakQueueItem(timestamps, newdata);
            lock (((ICollection)dataWritingQueue).SyncRoot)
            {
                dataWritingQueue.Enqueue(quI);
                ql = dataWritingQueue.Count;
            }
            dataWritingSyncEvent.NewItemEvent.Set();
            if (displayCount == 10)
            {
                setTextBox2(data.QueueLen.ToString());
                setTextBox3(ql.ToString());
            }
            if (gquitting)
            {
                return;
            }

        }



        //===========================GUI=====================

        public void setTextBox1(string o)
        {
            if (textBox1.InvokeRequired)
            {
                textBox1.BeginInvoke(
                   new MethodInvoker(
                   delegate () { setTextBox1(o); }));
            }
            else
            {
                textBox1.Text = o;
            }

        }

        public void setTextBox2(string o)
        {
            if (textBox2.InvokeRequired)
            {
                textBox2.BeginInvoke(
                   new MethodInvoker(
                   delegate () { setTextBox2(o); }));
            }
            else
            {
                textBox2.Text = o;
            }
        }

        public void setTextBox3(string o)
        {
            if (textBox3.InvokeRequired)
            {
                textBox3.BeginInvoke(
                   new MethodInvoker(
                   delegate () { setTextBox3(o); }));
            }
            else
            {
                textBox3.Text = o;
            }
        }
        public void setTextBox4(string o)
        {
            if (textBox4.InvokeRequired)
            {
                textBox4.BeginInvoke(
                   new MethodInvoker(
                   delegate () { setTextBox4(o); }));
            }
            else
            {
                textBox4.Text = o;
            }
        }


        /*
         *  Sets the size of the box and the corresponding locations of all relevant graphical components
         */
        public void SetSize()
        {
            int iS, lay, coy;
            iS = 10;
            lay = 4;
            coy = 20;
            label1.Location = new Point(iS, lay);
            textBox1.Location = new Point(iS, coy); iS = iS + textBox1.Width;
            label2.Location = new Point(iS, lay);
            textBox2.Location = new Point(iS, coy); iS = iS + textBox2.Width;
            label3.Location = new Point(iS, lay);
            textBox3.Location = new Point(iS, coy); iS = iS + textBox3.Width;
            label4.Location = new Point(iS, lay);
            textBox4.Location = new Point(iS, coy); iS = iS + textBox4.Width;
            label5.Location = new Point(iS, lay);
            numericUpDown1.Location = new Point(iS, coy); iS = iS + numericUpDown1.Width;
            label6.Location = new Point(iS, lay);
            numericUpDown2.Location = new Point(iS, coy); iS = iS + numericUpDown2.Width;
            label7.Location = new Point(iS, lay);
            numericUpDown3.Location = new Point(iS, coy); iS = iS + numericUpDown3.Width;
            buRecord.Location = new Point(ClientRectangle.Width - buRecord.Size.Width - 10, 20);
            buGraph.Location = new Point(ClientRectangle.Width - buGraph.Size.Width - 10, ClientRectangle.Height - buGraph.Size.Height - 20);

        }


        private void Form1_FormClosing(object sender, FormClosingEventArgs e)
        {
            int i;
            gquitting = true;
            cameraThread.Abort();
            myCamera.stopFrameGrab(cameraType);
            dataWritingSyncEvent.ExitThreadEvent.Set();
            for (i = 0; i < consumerd.Length; i++)
            {
                consumerd[i].mySyncEvent.ExitThreadEvent.Set();
            }

            for (i = 0; i < consumerd.Length; i++)
            {
                consumerd[i].myThread.Join();
            }
        }

        private void Form1_Resize(object sender, EventArgs e)
        {
            SetSize();
        }
        private void buGraph_Click(object sender, EventArgs e)
        {
            if (Application.OpenForms.OfType<Form2>().Count() == 0)
            {
                graphWindow = new Form2();
                graphThread = new Thread(Program.Main2);
                graphThread.SetApartmentState(ApartmentState.STA);
                graphThread.Start();
            }
            else
            {
                graphWindow.BringToFront();
            }
        }
        private void buRecord_Click(object sender, EventArgs e)
        {
            if (gRecord == false)
            {
                Frameco = 0;
                DTFrameCo0 = DateTime.Now;
                dayFrameCo0 = DayNr.GetDayNr(DTFrameCo0);


                buRecord.Text = "Stop Rec.";
                myDataWriter = new DataWriting();
                gRecord = true;
            }
            else
            {

                buRecord.Text = "Record";
                gRecord = false;
                myDataWriter.stopit();
            }
        }

        private void numericUpDown3_ValueChanged(object sender, EventArgs e)
        {
            gNoPeaks = (int)numericUpDown3.Value;
        }

    }
}
//}