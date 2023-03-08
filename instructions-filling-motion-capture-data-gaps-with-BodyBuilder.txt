%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%									     %%%%
%%%% Using Bodybuilder (BB) to fill the data gaps at start and end of trials %%%%
%%%%									     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Date written 15/09/17 by Franz Tapia Chaca

Works for standing trials (CPRT confirmed, BandB assumed)

1. Download BodyBuilder through the Vicon website
--------------------------------------------------

Link: https://www.vicon.com/downloads/software/bodybuilder

Version 3.6.4 is compatible with the dongle for Nexus 1.8.5, and if the website updates the version, then possibly try their archive (https://www.vicon.com/downloads/software/bodybuilder)

2. Learning to use BodyBuilder environment
-------------------------------------------

BodyBuilder and Nexus have similar environments and their uses overlap too (e.g. filling gaps, labelling, etc.). Though because it is different, the manual helps in understanding it. Uploaded on this folder 'Bodybuilder codes' is the mots recent manual found, but Vicon offers an older manual by 5 years on their website: https://www.vicon.com/downloads/documentation/bodybuilder-for-biomechanics

To understand completely how to fill data gaps at start and ends of trials (and be comfortable with labelling markers on BB when necessary (see Note 1)), it is suggested to read pages:
* 12 - 31
* 75 - 76 ('Using Kinetic Data: Terms and Definitions' not needed)
* 78 - 79 (from 'Applying the Model')
* 80 - 90 (stop before 'Antiflip' - I didn't read this as I didn't see a need)
* 110 - 119 (a similar problem to ours is explained in page 119)
* 123 - 139 (from 'Model script' - can come in helpful to learn syntax)

Note 1: When Nexus offers problems with labelling (i.e. markers are not continuous through frames, so one has to manually label every marker of every frame), BB surpasses this as the markers on the BB environment are continuous.

3. Writing the Model Scripts (.mod)
------------------------------------

Use Notepad or Notepad ++ (latter looks better)

There are different ways of going about this - my method uses the static standing trial data for referencing of the location of the missing markers.

Example: filling EPICEM gaps at start and end frames

++++++++++++
++ Code 1 ++
++++++++++++

{* Segment definition of Arm with LEPICEL as the Origin *}
Arm = [LEPICEL,LARM1-LEPICEL,LARM2-LEPICEL,yzx]			(line 1)

{* Local LEPICEM location from LEPICEL *}
{* Calculating average local location over all frames *}
{*   with respect to Arm segment *}
%LEPICEMavgArm=AVERAGE(LEPICEM/Arm)				(line 2)
PARAM(%LEPICEMavgArm)						(line 3)


Explanation:
This code is meant to run on the Static Standing Trial.
Line 1 creates a segment (cluster) of markers with its own coordinate system. LEPICEL is the origin and the y-axis falls along the line from LEPICEL-to-LARM1. The z-axis is made perpendicular to both lines LARM1-LEPICEL and LARM2-LEPICEL, and the x-axis is then calculated accordingly. The name of this segment is Arm.

Line 2 finds the location of LEPICEM in the coordinate frame of ARM in a number of steps:
1) LEPICEM/Arm: Convert the global (Nexus coordinate frame) coordinates of LEPICEM into local (Arm coordinate frame) coordinates with respect to the origin of Arm, using '/'.
2) AVERAGE(LEPICEM/Arm): Calculate the average location of the LEPICEM with respect to Arm. Function carried out over all frames.
3) %LEPICEMavgArm = AVERAGE(LEPICEM/Arm): Save this number on the variable with name %LEPICEMavgArm, wherein the '%' represents 'local coordinates'

Line 3 writes this variable onto the Parameter File (.mp) which every session folder has, and is accessed by BB every time a .mod script is run.


++++++++++++
++ Code 2 ++
++++++++++++

IF ((SAMPLE==FIRSTSAMPLE) OR (SAMPLE==LASTSAMPLE) OR (SAMPLE==(LASTSAMPLE-1))) THEN		(line 1)
	IF EXIST(LEPICEM)									(line 2)
		{* Do nothing *}
	ELSE											(line 3)
		{* Segment Definition of Arm *}
		Arm = [LEPICEL, LARM1-LEPICEL, LARM2-LEPICEL,yzx]				(line 4)
		%LEPICELorigin=(LEPICEL/Arm)							(line 5)

		{* Convert local marker to global marker, with respect to Arm segment *}
		LEPICEM=((%LEPICEMavgArm+%LEPICELorigin)*Arm)					(line 6)
		OUTPUT(LEPICEM)									(line 7)
	ENDIF											(line 8)
ELSE												(line 9)
	{* Do nothing *}
ENDIF												(line 10)

Explanation:
This code is meant to run on the dynamic standing trial that you want to gap fill for LEPICEM.
Line 1, 9, 10: SAMPLE is the current frame that the script is being run on. So, if the current frame is NOT the first frame, the last frame or the penultimate frame, ignore. But if YES, then run the following if-statement.
Line 2, 3, 8: If the marker LEPICEM exists on the current frame (i.e. not missing), then do nothing. But if NON-EXISTENT, then run following code.
Line 4: Redefine the segment Arm (as we are running this code on a different trial)
Line 5: Find the local coordinates, with respect to Arm segment, of LEPICEL, which is the origin. It can be expected that (0,0,0) is returned, but this was found needed for the next step.
Line 6: Add the average local coordinate of LEPICEM to the local origin, and convert this a global coordinate, using '*'.
Line 7: Show the new coordinates as a virtual, unlabelled marker on the workspace.


4. Running Model Scripts (.mod)
--------------------------------

Code 1:
- For the session that you want to fill a data gap, open the respective Static Standing Trial on BB
- Run .mod script
- Save parameter (.mp) file

Code 2:
- For that same session, open the desired dynamic standing trial
- Run .mod script
- Rotate the workspace to verify correct execution
- Click 'Save' to overwrite .c3d file.
