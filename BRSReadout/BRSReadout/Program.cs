﻿/**     BRSReadout
 **     
 **     Michael Ross and Krishna Venkateswara
 **     Eot-Wash Group
 **     Center for Experimental Nuclear Physics and Astrophysics (CENPA)
 **     University of Washington
 **
 **     Developed Winter 2015
 **
 **     Based on: BRSReadout by Trevor Arp
 **
 **     This program is the Data Aquisition software for a fringe autocollimator to be part of a tiltWash sensor
 **     being developed for the LIGO project.
 **     
 **     Last Updated: 1/15/16
 **/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;
using System.Diagnostics;

namespace BRSReadout
{

    /*
    *  Program class, this class contain the main method. It initializes the GUI and then runs the program.
    */
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            System.Diagnostics.Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.High;
            Application.Run(new Form1());
        }
        public static void Main2()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            System.Diagnostics.Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.High;
            Application.Run(new Form2());
        }
    }
}