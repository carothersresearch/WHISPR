

############################################################

#Input file
infile = '/Users/ryancardiff/Documents/GitHub/WHISPR/Experiments/230220_tdcb/230224_KOH_all_traces.csv' #Should be the .csv export from mass hunter - no modifications necessary.

#Output files
outraw = '/Users/ryancardiff/Documents/GitHub/WHISPR/Experiments/230220_tdcb/output_raw.csv' #Optional output file (.csv) that converts the mass hunter export into X Y1 Y2 ... Yn CSV
outmod = '/Users/ryancardiff/Documents/GitHub/WHISPR/Experiments/230220_tdcb/output_mod.csv' #Optional output file (.csv) of the subset of data that is being plotted (with scaling)
outwat = '/Users/ryancardiff/Documents/GitHub/WHISPR/Experiments/230220_tdcb/output.pdf' #Output figure location, should be .pdf to import into illustrater if needed.

############################################################
#Console output
print_sname = True #Print sample name to the console for easy identification

############################################################

#PLOT PARAMETERS
y_reorder = [] #Note: Comma required between numbers!
#What this does: Re order columns or chose only a subset of the data to plot. If empty ( = []), use all columns in given order. Index starts at 0. 
#Example: [0 , 2, 3 , 9 , 5] will only plot the y0, y2, y3, y9, y5 columns in that order. Use print_sname = True to see sample labels. Note: Comma required between numbers!

#Title
title_text = '' #Optional - Leave empty if no title desired.
title_size = 20 #Font size of title text

#Y-Axiy_padding = 1.5 #Scaling factor for chromatograms (e.g. 1.5; larger = shrink). Used to make things look neater.
y_padding = 1.5

#X-Axis
t_min = 3 #Min time in units data is stored in (seconds, minutes, hours)
t_max = 5 #Max time in units data is stored in (seconds, minutes, hours)
x_min = 2 #Minimum on x-axis
x_max = 5 #Maximum on x-axis
x_tick = 0.5 #Incremenet for x-tick marks
x_tick_size = 16 #Font size for x-tick marks
x_label = 'Time (min)' #X-axis title. Leave blank ( = '') if no title desired.
x_label_size = 18 #Font size for X-title

#Lines
axis_linewidth = 1.5 #linewidth of the X and Y axis
chrom_linewidth = 1.5 #linewidth of the chromatograms

############################################################

