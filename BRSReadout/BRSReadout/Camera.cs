﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Windows.Forms;
using System.Runtime.InteropServices;
using System.Diagnostics;
using System.IO;
using System.Threading;
using PylonC.NET;
using TwinCAT.Ads;
using BRSReadout;


/*
 *  Camera Class aquires data coming for the Basler line camera. Must be run in seperate thread to allow for continuous image acquisition.
 */
class Camera
{
    public delegate void frameCallbackDelegate(PylonBuffer<Byte> buffer);
    public delegate void dataDelegate(ushort[] data);

    dataDelegate masterDataDelegate;
    public ushort[] data = new ushort[4096];

    PYLON_DEVICE_HANDLE hDev = new PYLON_DEVICE_HANDLE();

    uint numDevices;
    PYLON_STREAMGRABBER_HANDLE hGrabber;
    PYLON_WAITOBJECT_HANDLE hWait;
    uint payloadSize;
    Dictionary<PYLON_STREAMBUFFER_HANDLE, PylonBuffer<Byte>> buffers;
    PylonGrabResult_t grabResult;

    uint nStreams;
    bool isAvail;
    bool isReady;

    PYLON_DEVICE_INFO_HANDLE prop;
    string ip;
    uint j;



    public Camera()
    {
    }
    public void cameraInit()
    {
        //Initiates camera to run in 12 bit mode with continuous acquisition.
        int i = 0;
        uint NUM_BUFFERS = 10;
        Pylon.Initialize();
        numDevices = Pylon.EnumerateDevices();
        while (ip != "10.110.136.138")
        {
            hDev = Pylon.CreateDeviceByIndex(j);
            prop = Pylon.DeviceGetDeviceInfoHandle(hDev);
            ip = Pylon.DeviceInfoGetPropertyValueByIndex(prop, 8);
            j++;
            if (j > 100) { break; }
        }
        try
        {
            Pylon.DeviceOpen(hDev, Pylon.cPylonAccessModeControl | Pylon.cPylonAccessModeStream);
        }
        catch (Exception ex)
        {
            EmailError.emailAlert(ex);
            throw (ex);
        }

        isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_PixelFormat_Mono12");
        Pylon.DeviceFeatureFromString(hDev, "PixelFormat", "Mono12");
        Pylon.DeviceFeatureFromString(hDev, "AcquisitionMode", "Continuous");
        Pylon.DeviceSetIntegerFeature(hDev, "Height", 1);
        Pylon.DeviceSetFloatFeature(hDev, "ExposureTimeAbs", 500); //Exposure time is in microseconds and rounded to the closest 100 ns.
        isAvail = Pylon.DeviceFeatureIsWritable(hDev, "GevSCPSPacketSize");
        if (isAvail)
        {
            Pylon.DeviceSetIntegerFeature(hDev, "GevSCPSPacketSize", 1500);
        }
        payloadSize = checked((uint)Pylon.DeviceGetIntegerFeature(hDev, "PayloadSize"));
        nStreams = Pylon.DeviceGetNumStreamGrabberChannels(hDev);
        hGrabber = Pylon.DeviceGetStreamGrabber(hDev, 0);
        Pylon.StreamGrabberOpen(hGrabber);
        hWait = Pylon.StreamGrabberGetWaitObject(hGrabber);
        Pylon.StreamGrabberSetMaxNumBuffer(hGrabber, NUM_BUFFERS);
        Pylon.StreamGrabberSetMaxBufferSize(hGrabber, payloadSize);
        Pylon.StreamGrabberPrepareGrab(hGrabber);
        buffers = new Dictionary<PYLON_STREAMBUFFER_HANDLE, PylonBuffer<Byte>>();
        for (i = 0; i < NUM_BUFFERS; ++i)
        {
            PylonBuffer<Byte> buffer = new PylonBuffer<byte>(payloadSize, true);
            PYLON_STREAMBUFFER_HANDLE handle = Pylon.StreamGrabberRegisterBuffer(hGrabber, ref buffer);
            buffers.Add(handle, buffer);
        }
        i = 0;
        foreach (KeyValuePair<PYLON_STREAMBUFFER_HANDLE, PylonBuffer<Byte>> pair in buffers)
        {
            Pylon.StreamGrabberQueueBuffer(hGrabber, pair.Key, i++);
        }
        Pylon.DeviceExecuteCommandFeature(hDev, "AcquisitionStart");
    }

    public void startFrameGrab(int nr, int trigmode, dataDelegate dd)
    {
        //Start frame grabbing loop and pushes data into the data delegate as a ushort array.
        int bufferIndex;
        while (true)
        {
            try
            {
                PylonBuffer<Byte> buffer;
                Pylon.WaitObjectWait(hWait, 100);

                Pylon.StreamGrabberRetrieveResult(hGrabber, out grabResult);

                bufferIndex = (int)grabResult.Context;
                while (buffers.TryGetValue(grabResult.hBuffer, out buffer) != true) ;
                masterDataDelegate = dd;
                Buffer.BlockCopy(buffer.Array, 0, data, 0, buffer.Array.Length);
                masterDataDelegate(data);
                Pylon.StreamGrabberQueueBuffer(hGrabber, grabResult.hBuffer, bufferIndex);
            }
            catch (Exception ex)
            {
                EmailError.emailAlert(ex);
                throw (ex);
            }
        }
    }

    public void stopFrameGrab()
    {
        //Stops acquisition and closes camera.
        Pylon.DeviceExecuteCommandFeature(hDev, "AcquisitionStop");

        Pylon.StreamGrabberCancelGrab(hGrabber);

        do
        {
            isReady = Pylon.StreamGrabberRetrieveResult(hGrabber, out grabResult);

        } while (isReady);

        foreach (KeyValuePair<PYLON_STREAMBUFFER_HANDLE, PylonBuffer<Byte>> pair in buffers)
        {
            Pylon.StreamGrabberDeregisterBuffer(hGrabber, pair.Key);
            pair.Value.Dispose();
        }
        buffers = null;

        Pylon.StreamGrabberFinishGrab(hGrabber);

        Pylon.StreamGrabberClose(hGrabber);

        Pylon.DeviceClose(hDev);
        Pylon.DestroyDevice(hDev);

        Console.ReadLine();

        Pylon.Terminate();

    }

}



