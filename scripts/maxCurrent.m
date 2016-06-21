rtn = vicon_process_trial('~/research/data/PA_A01_FF_Jan_22_2015/PA_A02_SS_F1_02')


t1 = min(rtn.id.RAnkleMoment.X)
t2 = min(rtn.id.LAnkleMoment.X)

currentAtPeakDifference = (t1 - t2)*85*0.4
