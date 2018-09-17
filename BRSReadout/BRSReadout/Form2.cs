using LiveCharts;
using LiveCharts.Wpf;
using System;
using System.Collections.Generic;
using System.Configuration;
using System.Drawing;
using System.Threading;
using System.Windows.Forms;

namespace BRSReadout
{
    public partial class Form2 : Form
    {
        public static object graphLock = Form1.graphLock;
        public static Queue<Form1.graphData> graphQueue = Form1.graphQueue;
        public static AutoResetEvent graphSignal = Form1.graphSignal;
        public bool run=true;
        Thread dataLoopThread;
        public void stopLoop()
        {
            run = false;
        }
        public void dataLoop()
        {
            while (true)
            {
                graphSignal.WaitOne();
                Form1.graphData inData = new Form1.graphData();
                do
                {
                    inData.frame = null;
                    lock (graphLock)
                    {
                        if (graphQueue.Count > 0)
                        {
                            inData = graphQueue.Dequeue();
                        }
                        if (inData.frame != null)
                        {
                            updatePlot(inData.frame, inData.angle);
                        }
                    }
                }
                while (inData.frame != null);
            }
        }
        public Form2()
        {
            InitializeComponent();
            ushort[] rawFrame = Form1.frame;
            double[,] data = Form1.newdata;
            imagePlot.Series = new SeriesCollection
            {
                new LineSeries
                {
                    Title = "Image",
                    Values = new ChartValues<int> (),
                    Fill=System.Windows.Media.Brushes.Transparent,
                    Stroke = System.Windows.Media.Brushes.Red,
                    LineSmoothness = 0,
                    PointGeometry = null,
                    PointGeometrySize = 0
                }
            };
            anglePlot.Series = new SeriesCollection
            {
                new LineSeries
                {
                    Title = "Angle",
                    Values = new ChartValues<double>(),
                    Fill=System.Windows.Media.Brushes.Transparent,
                    LineSmoothness = 0,
                    PointGeometry = null,
                    PointGeometrySize = 0
                }
            };

            imagePlot.AxisX.Add(
            new Axis
            {
                MinValue = 0,
                MaxValue = 4096
            });

            imagePlot.AxisY.Add(
            new Axis
            {
                MinValue = 0,
                MaxValue = Math.Pow(2, double.Parse(ConfigurationManager.AppSettings.Get("cameraBitDepth")))
            });


            imagePlot.DisableAnimations = true;
            imagePlot.Hoverable = false;
            imagePlot.DataTooltip = null;

            anglePlot.DisableAnimations = true;
            anglePlot.Hoverable = false;
            anglePlot.DataTooltip = null;

            imagePlot.Anchor = AnchorStyles.Left;
            anglePlot.Anchor = AnchorStyles.Right;
            imagePlot.Anchor = AnchorStyles.Top;
            anglePlot.Anchor = AnchorStyles.Top;

            if (dataLoopThread == null)
            {
                dataLoopThread = new Thread(dataLoop);
                dataLoopThread.Start();
            }
            else
            {
                dataLoopThread.Start();
            }

        }
        public static IEnumerable<T> ToEnumerable<T>(Array target)
        {
            foreach (var item in target)
            {
                yield return (T)item;
            }
        }
        public void updatePlot(ushort[] rawFrame, double data)
        {
            if (anglePlot.Series[0].Values.Count > 1000)
            {
                anglePlot.Series[0].Values.RemoveAt(0);
                anglePlot.Series[0].Values.Add(data);
            }
            else
            {
                anglePlot.Series[0].Values.Add(data);
            }

            imagePlot.Series[0].Values.Clear();
            var inFrame = Array.ConvertAll(rawFrame, item => (int)item);
            imagePlot.Series[0].Values.AddRange(ToEnumerable<object>(inFrame));
        }
        public void showWindow()
        {
            this.Show();
            this.BringToFront();
            try
            {
                dataLoopThread.Start();
            }
            catch { }
        }
        private void Form2_FormClosing(object sender, FormClosingEventArgs e)
        {
            stopLoop();
            this.Hide();
            e.Cancel=true;
        }

    }
}
