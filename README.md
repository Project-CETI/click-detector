# detector
"# Coda_detector" 

In order for the Coda detector to work, please:
1) Download the toolbox: "MVTB-4.3" in the following link: https://petercorke.com/toolboxes/machine-vision-toolbox/
2) From within the MATLAB file browser double click on this file, it will install and configure the paths correctly

---

The script applies two detection schemes, one is customized for recordings from tags, and the other for recordings from the buoy.
The control over the detection choice is made in line 12 via: Tag_flag (set 1 for Dtag, 0 for Buoy).

For Buoy data: 
The contorl over the click type detector is made in line 10 via: Detector_flag (1- apply coda detector | 0- apply echolocation clicks detector)

To visualize the detector outputs: set  Plot_flag=1 (line 11)

The time of arrivals of the detected click will be automatically saved in an .xls file with the same name as that of the recording.
