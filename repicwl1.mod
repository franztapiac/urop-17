{* Gap filling of REPICWL: Part 1 *}
{* Reading from static standing trial *}

CPalm = ((REPICWM+RPALMM+RPALML+RPALMC)/4)
RightPalm = [CPalm, RPALMC-REPICWM, RPALML-REPICWM,yzx]

%REPICWLavgPalm=AVERAGE(REPICWL/RightPalm)
PARAM(%REPICWLavgPalm)

