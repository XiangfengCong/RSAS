

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                    */
/*                                              Version RCS_REG V1.5 beta                                             */
/*                                              -------------------------                                             */
/*                                                                                                                    */
/*                                                     CONTACT INFO                                                   */
/*                                                     ------------                                                   */
/*                                                                                                                    */
/*  Loic Desquilbet, PhD                                                                                              */
/*  E-mail: loic.desquilbet@gmail.com                                                                                 */
/*                                                                                                                    */
/*                                                                                                                    */
/*                                                     INTRODUCTION                                                   */
/*                                                     ------------                                                   */
/*                                                                                                                    */
/* The RCS_Reg SAS macro performs two major tasks: (1) it creates a restrictive cubic spline (RCS) function of a      */
/* continuous exposure with 3-5 knots; (2) it displays the dose-response association (with its 95% confidence         */
/* interval) between the continuous exposure and an outcome, using linear, logistic, or Cox regressions, as well as   */
/* GEE linear and logistic models. To display the dose-response association, the macro uses a reference value         */
/* (assigned by the macro or chosen by the user); it displays estimated differences (in case of a linear              */
/* regression), Log(Odds Ratio) (in case of a logistic regression), or Log(Hazard Ratio) (in case of a Cox model),    */
/* when comparing subjects with any value of the continuous exposure with subjects with the reference value of the    */
/* continuous exposure.                                                                                               */
/*                                                                                                                    */
/*                                                MEANING OF THE PARAMETERS                                           */
/*                                                -------------------------                                           */
/*                                                                                                                    */
/*                                                                                                                    */
/* The "**" sign indicates that the parameter must always be specified. The "*" sign indicates that the parameter     */
/* must be specified if a regression model is requested. Otherwise, the parameter is optional. In this case and if    */
/* the parameter is not specified, the default value is used.                                                         */
/*                                                                                                                    */
/* INFILE**: name of the original SAS datafile that contains the continuous exposure(s) that one wants to code        */
/* using a restricted cubic spline (RCS) function. A libname must NOT precede the name of this datafile. If the       */
/* original SAS datafile is located in a specific file, use the DIR_DATA parameter (see below).                       */
/*                                                                                                                    */
/* OUTFILE: name of the new SAS datafile, created in the SAS Working library, that will contain the splines of the    */ 
/* continuous exposure(s). If not specified, the name of this outfile datafile is the name of the original SAS        */
/* datafile followed by "_rcs".                                                                                       */
/*                                                                                                                    */
/* MAIN_SPLINE_VAR**: name of the main continuous exposure that one wants to code using a RCS function. The main      */
/* continuous exposure is the continuous exposure for which one (optionally) wants to display the curve of the        */
/* dose-response association with the outcome. If the name of this variable contains more than 12 characters, its     */
/* name will be truncated to 11 characters, and the 12th character will be "_".                                       */
/*                                                                                                                    */
/* AVK_MSV: by default, the knots of the main continuous exposure are defined according to the percentiles of its     */
/* distribution (AVK_MSV=0). If AVK_MSV is set to 1, the user chooses to define Arbitrary Values of the Knots.        */
/*                                                                                                                    */
/* KNOTS_MSV: if AVK_MSV=0, KNOTS_MSV is the list of the values of the percentiles of the distribution of the         */
/* main continuous exposure. These values must be separated by a space and be provided in the ascendant order         */
/* (the lowest percentile comes first). For instance, if AVK_MSV=0, "KNOTS_MSV=5 25 75 95" means that the user        */
/* chooses the values of the 5th, 25th, 75th, and 95th percentiles of the main continuous exposure as the values      */
/* of the knots. If AVK_MSV=0, the default list is "5 50 95". If AVK_MSV=1, KNOTS_MSV must be specified, and it is    */
/* the list of arbitrary values of the main continuous exposure that define the knots. For instance, if the main      */
/* continuous exposure is the age of individuals, and if AVK_MSV=1, "KNOTS_MSV=23 29 45" means that the knots are     */
/* located at 23, 29, and 45 years old.                                                                               */
/*                                                                                                                    */
/* OTH_SPLINE_VAR1: first other continuous exposure that one wants to include into the regression model using a       */
/* RCS function, as a potential confounder. If the name of this variable contains more than 12 characters, its name   */
/* will be truncated to 11 characters, and the 12th character will be "_".                                            */
/*                                                                                                                    */
/* AVK_OSV1: same as AVK_MSV for the first other continuous exposure that is specified in OTH_SPLINE_VAR1 parameter.  */
/*                                                                                                                    */
/* KNOTS_OSV1: same as KNOTS_MSV for the first other continuous exposure that is specified in OTH_SPLINE_VAR1         */
/* parameter.                                                                                                         */
/*                                                                                                                    */
/* OTH_SPLINE_VAR2-OTH_SPLINE_VAR10: second to (max) tenth other continuous exposures that one wants to include       */
/* into the regression model using a RCS function, as potential confounders.                                          */
/*                                                                                                                    */
/* KNOTS_OSV2-KNOTS_OSV10: same as KNOTS_OSV1 for the second to the (max) tenth other continuous exposures that are   */
/* specified in OTH_SPLINE_VAR2-OTH_SPLINE_VAR10 parameters.                                                          */
/*                                                                                                                    */
/* AVK_OSV2-AVK_OSV10: same as AVK_OSV1 for the second to the (max) tenth other continuous exposures that are         */
/* specified in OTH_SPLINE_VAR2-OTH_SPLINE_VAR10 parameters.                                                          */
/*                                                                                                                    */
/* DIR_DATA: name of the input directory that contains the original SAS datafile specified in the INFILE parameter.   */
/* The name must NOT be typed between quotes (for instance, DIR_DATA=C:\SAS\data). By default, the input directory is */
/* the SAS Working library.                                                                                           */
/*                                                                                                                    */
/* TYP_REG: specifies the type of the regression. Assign "lin" (without quotes) to perform a linear regression,       */
/* "log" for a logistic regression, or "cox" for a Cox model. If TYP_REG is not specified, the macro only creates     */
/* the splines of the continuous exposure(s) in the outfile datafile specified in the OUTFILE parameter, without      */
/* performing any regression.                                                                                         */
/*                                                                                                                    */
/* DEP_VAR*: name of the dependent variable (outcome) in the regression model. For a linear regression, DEP_VAR       */
/* must be continuous. For a logistic regression, DEP_VAR must value 1 for cases and 0 for controls. For a Cox model, */
/* DEP_VAR must value 1 for events and 0 for censored individuals. If the name of this variable contains more than    */
/* 16 characters, its name will be truncated to 15 characters, and the 16th character will be "_".                    */
/*                                                                                                                    */
/* SURV_TIME_VAR: name of the survival time variable in the case of a Cox model. SURV_TIME_VAR must be specified      */
/* if TYP_REG=cox is specified.                                                                                       */
/*                                                                                                                    */
/* ADJUST_VAR: list of the adjusted variables included into the regression model, except spline variables specified   */
/* in the OTH_SPLINE_VARk parameters (k   {1, …, 10}). The list elements must be separated by a space. If a           */
/* continuous exposure is specified in this ADJUST_VAR parameter, its association with the outcome is assumed to be   */
/* linear. If ADJUST_VAR and OTH_SPLINE_VARk are not specified, the regression will display an unadjusted             */
/* dose-response association between the main continuous exposure and the outcome.                                    */
/*                                                                                                                    */
/* PRGM_STATEMENTS_COX: in the case of a Cox model, the macro uses the PROC PHREG SAS procedure. This procedure       */
/* enables to create new explanatory variables or modify the values of explanatory variables in the MODEL statement   */
/* of the procedure. The parameter PRGM_STATEMENTS_COX can be used to do this/theses same task(s) when performing a   */
/* Cox model with the macro (for instance, to create time by covariate interactions).                                 */
/*                                                                                                                    */
/* BY_FACTOR: name of a categorical variable (which must be numerical) used for performing a "by factor" analysis     */
/* according to this variable, i.e., an analysis stratified on this categorical variable. If a "by factor" analysis   */
/* is requested, the name of the outfile datafile (specified in the OUTFILE parameter) will be followed by "_k",      */
/* where k is the rank of the category (the lowest category being the first); there will therefore be as many outfile */
/* datafiles as the number of categories of the categorical variable.                                                 */
/*                                                                                                                    */
/* MISSING_BF: if set to 1, a missing value on the variable specified in BY_FACTOR parameter will be treated as a     */
/* non missing value (and assigned to -99). In this case, the first outfile datafile will be for the missing values   */
/* of the categorical variable specified in the BY_FACTOR parameter. By default, the analysis is not performed for    */
/* missing data on this variable (MISSING_BF=0).                                                                      */
/*                                                                                                                    */
/* SUBJECT_VAR: in case of multiple records per subject, the macro can perform linear or logistic GEE models by       */
/* using the "REPEATED SUBJECT" command in PROC GENMOD SAS procedure. To perform a GEE model (and only to do so),     */
/* assign the name of the subject variable.                                                                           */
/*                                                                                                                    */
/* WORK_CORR_MATRIX: structure of the working correlation matrix in case of GEE models, that will be used by the      */
/* PROC GENMOD procedure. A list of available matrix structures can be found at                                       */
/* http://support.sas.com/onlinedoc/913/docMainpage.jsp. For instance, WORK_CORR_MATRIX=mdep(3) requests a            */
/* 3-dependent working correlation matrix. This parameter must be specified if SUBJECT_VAR parameter is specified.    */
/*                                                                                                                    */
/* REF_VAL: reference value of the main continuous exposure that will be used to display the dose-response            */
/* association with the outcome when comparing individuals with any value of the main continuous exposure with        */
/* individuals with the value assigned in REF_VAL. By default, REF_VAL is the median of the main continuous exposure. */  
/*                                                                                                                    */
/* SPECIF_VAL: the macro can provide in the SAS OUTPUT window the estimated value of the association (with its        */
/* 95% confidence interval) between the main continuous exposure and the outcome for one or more specific values of   */
/* the exposure (assigned as a list in SPECIF_VAL) when compared with the reference value (assigned in REF_VAL).      */
/* The list of these specific values must be separated by a space. If SPECIF_VAL is specified, a dataset named        */
/* "List_specif_values" is created in the SAS Working library. If a "by factor" analysis is requested, the name of    */
/* this datafile will be followed by "_k", where k is the rank of the category of the variable specified in           */
/* BY_FACTOR parameter (the lowest category being the first).                                                         */
/*                                                                                                                    */
/* ROUND: specifies if the value of the association provided by using SPECIF_VAL parameter has to be rounded. For     */
/* instance, if set to 0.01, it means that the value of the association (and the ones in the 95% confidence interval) */
/* will be rounded to the nearest 0.01 value when displayed in the SAS OUTPUT window. By default, values are not      */
/* rounded.                                                                                                           */
/*                                                                                                                    */
/* EXP_BETA: when set to 0, the macro displays Ln(ORs) or Ln(HRs); when set to 1, the macro displays ORs or HRs.      */
/* The default value is 0. Of note, this parameter has no effect if TYP_REG=lin (linear regression).                  */
/*                                                                                                                    */
/* PRINT_OR_HR: if set to 1 and if EXP_BETA is set to 1, the macro provides in the SAS OUTPUT window ORs or HRs       */
/* (with their 95% confidence intervals). In this case, a datafile named "List_or_hr" is created in the SAS           */
/* Working library. If a "by factor" analysis is requested, the name of this datafile will be followed by "_k",       */
/* where k is the rank of the category of the variable specified in BY_FACTOR parameter (the lowest category being    */
/* the first). By default, ORs or HRs are not provided (PRINT_OR_HR=0). This parameter has no effect if EXP_BETA=0.   */
/*                                                                                                                    */
/* WHERE: this parameter enables to create the splines and to perform regressions in a selection of individuals,      */
/* like the usual WHERE statement in all SAS procedures.                                                              */
/*                                                                                                                    */
/* HISTOGRAM: if set to 1, the macro displays the distribution of the main continuous exposure using an histogram     */
/* (HISTOGRAM statement in the PROC UNIVARIATE procedure). The default value is 0 (no histogram is displayed).        */
/*                                                                                                                    */
/* NO_GRAPH: if set to 1, no graph is displayed in the SAS OUTPUT window. The graph is displayed by default           */
/* (NO_GRAPH=0).                                                                                                      */
/*                                                                                                                    */
/* DISPLAY_KNOTS: if set to 1, the macro displays the knots on the dose-response curve as dots. If set to 0, the      */
/* knots are not displayed on the curve. The default value is 1.                                                      */
/*                                                                                                                    */
/* X_REF_LINE: if set to 1, a vertical dashed green line is displayed to materialize the reference value of the       */
/* main continuous exposure. There is no vertical line by default (X_REF_LINE=0).                                     */
/*                                                                                                                    */
/* Y_REF_LINE: if set to 1, a horizontal dashed green line is displayed to materialize the null hypothesis H0.        */
/* The Y-coordinate of the line is 0 when a linear regression is requested, or when a logistic or a Cox model         */
/* with EXP_BETA=0 is requested. The Y-coordinate of the line is 1 when a logistic or a Cox model with EXP_BETA=1     */
/* is requested. There is no horizontal line by default (Y_REF_LINE=0).                                               */
/*                                                                                                                    */
/* MIN_XAXIS: by default, the graph displaying the dose-response curve starts, on the X-axis, at the minimum value    */
/* of the main continuous exposure. Assigning a value to MIN_XAXIS enables the graph to start at this assigned value. */ 
/*                                                                                                                    */
/* MAX_XAXIS: by default, the graph displaying the dose-response curve ends, on the X-axis, at the maximum value of   */
/* the main continuous exposure. Assigning a value to MAX_XAXIS enables the graph to end at this assigned value.      */
/*                                                                                                                    */
/* NO_TITLE: if set to 1, there is no title on the graph. By default, there is a title (NO_TITLE=0).                  */
/*                                                                                                                    */
/* NO_LABEL_X: if set to 1, there is no label on the X-axis. By default, there is a label of the X-axis               */
/* (NO_LABEL_X=0), and it is the name of the main continuous exposure.                                                */
/*                                                                                                                    */
/* NO_LABEL_Y: if set to 1, there is no label on the Y-axis. By default, there is a label of the Y-axis               */
/* (NO_LABEL_Y=0), of which the text depends on the type of the regression.                                           */
/*                                                                                                                    */
/* NO_LEGEND: if set to 1, there is no legend under the X-axis. By default, the legend is below the X-axis            */
/* (NO_LEGEND=0).                                                                                                     */
/*                                                                                                                    */
/* PRINT_COVAR_MAT: if set to 1, the estimated regression parameter covariance matrix is displayed in the SAS         */ 
/* OUTPUT window. By default, this matrix is not displayed (PRINT_COVAR_MAT=0).                                       */
/*                                                                                                                    */
/* NO_DELETE_FILES: by default, the macro deletes all temporary files that have been created in the SAS Working       */
/* library (NO_DELETE_FILES=0). To keep these files, NO_DELETE_FILES must be assigned to 1.                           */
/*                                                                                                                    */
/*                                                                                                                    */
/*                                  MODIFICATIONS SINCE V1.0 (published in Stat Med, 2010)                            */
/*                                  ------------------------------------------------------                            */
/*                                                                                                                    */
/* The OUTPUT_PRED parameter has been removed, and the SAS macro now systematically provides the predicted value of   */
/* the outcome for each observation, in the outfile datafile (specified in OUTFILE parameter). For the linear         */
/* regression, the macro provides the predicted value of the continuous outcome as well as its lower and upper 95%    */
/* confidence limits; for the logistic regression, the macro provides the predicted value of the probability of the   */
/* binary outcome as well as its lower and upper 95% confidence limits; for the Cox model, the macro provides the     */
/* predicted value of the survival function S(t) computed by using the product-limit method (default method in the    */
/* PHREG SAS procedure).                                                                                              */
/*                                                                                                                    */
/* A bug on AIC display in the LOG window has been corrected.                                                         */
/*                                                                                                                    */
/* A bug with the WHERE statement has been corrected.                                                                 */
/*                                                                                                                    */
/* Now, the parameter DIR_DATA must NOT be typed between quotes (for instance, DIR_DATA=C:\SAS\data).                 */
/*																													  */
/* From v1.5, the macro can be run on SAS Studio Online (SAS University Edition) since PROC GPLOT has been replaced   */
/* by PROC SGPLOT. Therefore, some graphical characteristics may be different between v1.4 and v1.5. Please send me   */
/* an email if you have any troubles with plots. Of note, if you plan to use this macro on SAS Studio Online, you     */
/* must set-up the language of your Internet Browser in English.                                                      */
/*                                                                                                                    */
/*                                                                                                                    */
/*                                                                                                                    */
/*                                           ADDITIONAL PARAMETERS SINCE V1.0                                         */
/*                                           --------------------------------                                         */
/*                                                                                                                    */
/* PRINT_OBS: if set to 1, a "+" sign is displayed at the bottom of the graph for each observed value of the variable */
/* assigned in the MAIN_SPLINE_VAR parameter that has been used in the model (i.e., observations with missing data on */
/* variables included into the model are not displayed). By default, no "+" sign is displayed (PRINT_OBS=0).          */
/*                                                                                                                    */
/* WIDTH_BAND_OBS: width of the band of the "+" signs displaying observed values of the variable assigned in the      */
/* MAIN_SPLINE_VAR parameter, when PRINT_OBS=1. It is recommended that the value assigned in WIDTH_BAND_OBS should    */
/* not be greater than 5, although any positive value is allowed. The default value is 3.                             */
/*                                                                                                                    */
/* EXPORT_CURVES : If one wants to display the dose-response association using another sofware such as Excel,         */
/* specifying EXPORT_CURVES=1 exports a txt file (separator=TAB; name, "export_curves_rcs.txt"). If PRINT_OBS=1, two  */
/* columns corresponding the X and Y coordinates of the "+" signs are added. In this case, the first 201 lines must   */
/* be used to display the dose-response association with its confidence intervals, and the reminders must be used to  */
/* display the "+" signs, if the user wants so. EXPORT_CURVES=0 by default (no export of a txt file).                 */
/*                                                                                                                    */
/* DIR_EXPORT_CURVES : name of the export directory where one wants to export the txt file that will be used to       */
/* display the dose-response relationship with another software. This parameter must NOT be type between quotes (like */
/* DIR_DATA) and must be specified if EXPORT_CURVES=1. This parameter has no effect if EXPORT_CURVES=0 or if          */
/* EXPORT_CURVES is not specified.                                                                                    */
/*                                                                                                                    */
/*                                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------------ */


%macro RCS_Reg(	infile= ,
		outfile= ,
		main_spline_var= ,
		AVK_MSV= ,
		knots_MSV= ,
		oth_spline_var1= ,oth_spline_var2= ,oth_spline_var3= ,oth_spline_var4= ,oth_spline_var5= ,
		oth_spline_var6= ,oth_spline_var7= ,oth_spline_var8= ,oth_spline_var9= ,oth_spline_var10= ,
		AVK_OSV1= ,AVK_OSV2= ,AVK_OSV3= ,AVK_OSV4= ,AVK_OSV5= ,
		AVK_OSV6= ,AVK_OSV7= ,AVK_OSV8= ,AVK_OSV9= ,AVK_OSV10= ,
		knots_OSV1= ,knots_OSV2= ,knots_OSV3= ,knots_OSV4= ,knots_OSV5= ,
		knots_OSV6= ,knots_OSV7= ,knots_OSV8= ,knots_OSV9= ,knots_OSV10= ,
		dir_data= ,
		dep_var= ,
		surv_time_var= ,
		typ_reg= ,
		by_factor= ,
		missing_BF= ,
		adjust_var= ,
		prgm_statements_cox= ,
		subject_var= ,
		work_corr_matrix= ,
		ref_val= ,
		specif_val= ,
		round= , 
		exp_beta= ,
		print_OR_HR= ,
		where= ,
		histogram= ,
		no_graph= ,
		display_knots= ,
		X_ref_line= ,
		Y_ref_line= ,
		print_obs= ,
		width_band_obs= ,
		min_Xaxis= ,
		max_Xaxis= ,
		no_title= ,no_label_X= ,no_label_Y= ,no_legend= ,
		print_covar_mat= ,
		export_curves= ,
		dir_export_curves= ,
		no_delete_files= );


option nonotes;

%global AIC	FOR_RENAME 	I 	I_BIS 	INC	J 	K 	KNOT_INDEX 	KNOT1_MSV 	KNOT1_OSV1 	
		KNOT1_OSV2	KNOT1_OSV3	KNOT1_OSV4	KNOT1_OSV5	KNOT1_OSV6	KNOT1_OSV7	KNOT1_OSV8	KNOT1_OSV9	
		KNOT1_OSV10	KNOT2_MSV 	KNOT2_OSV1	KNOT2_OSV2	KNOT2_OSV3	KNOT2_OSV4	KNOT2_OSV5	KNOT2_OSV6	
		KNOT2_OSV7	KNOT2_OSV8	KNOT2_OSV9	KNOT2_OSV10	KNOT3_MSV 	KNOT3_OSV1 	KNOT3_OSV2	KNOT3_OSV3	KNOT3_OSV4	
		KNOT3_OSV5	KNOT3_OSV6	KNOT3_OSV7	KNOT3_OSV8	KNOT3_OSV9	KNOT3_OSV10	KNOT4_MSV 	KNOT4_OSV1	KNOT4_OSV2	
		KNOT4_OSV3	KNOT4_OSV4	KNOT4_OSV5	KNOT4_OSV6	KNOT4_OSV7	KNOT4_OSV8	KNOT4_OSV9	KNOT4_OSV10	KNOT5_MSV 	
		KNOT5_OSV1	KNOT5_OSV2	KNOT5_OSV3	KNOT5_OSV4	KNOT5_OSV5	KNOT5_OSV6	KNOT5_OSV7	KNOT5_OSV8	KNOT5_OSV9	
		KNOT5_OSV10	LEGEND 	LOG_L 	MAX_XAXIS_FA 	MEAN_MSV 	MEAN_OSV1 	MEAN_OSV2	MEAN_OSV3	
		MEAN_OSV4	MEAN_OSV5	MEAN_OSV6	MEAN_OSV7	MEAN_OSV8	MEAN_OSV9	MEAN_OSV10	MESSAGE_ADD1	
		MESSAGE_ADD2	MESSAGE_REF_VAL	MIN_XAXIS_FA 	MISTAKE 	NAME_MSV0 	NAME_MSV1 	NAME_MSV2 	NAME_MSV3 	
		NAME_OSV_CURR 	NAME10 	NAME11 	NB_CLASSES_BYFACT 	NB_KNOTS_CURR	NB_KNOTS_MINUS1 	NB_KNOTS_MSV 	
		NB_KNOTS_OSV1 	NB_KNOTS_OSV2	NB_KNOTS_OSV3	NB_KNOTS_OSV4	NB_KNOTS_OSV5	NB_KNOTS_OSV6	
		NB_KNOTS_OSV7	NB_KNOTS_OSV8	NB_KNOTS_OSV9	NB_KNOTS_OSV10	NB_OSV 	NB_PARAMETERS 	NB_SP_VALUES 	
		NO_CL_BYFACT 	OSV_CURR 	OUTFILE_FOR_ANALYSIS 	P_0_MSV	P_0_OSV1	P_100_MSV	P_100_OSV1	P_2_5_MSV 	
		P_2_5_OSV1 	P_20_MSV 	P_20_OSV1	P_25_MSV	P_25_OSV1	P_33_MSV 	P_33_OSV1	P_40_MSV	P_40_OSV1	
		P_5_MSV	P_5_OSV1 	P_50_MSV	P_50_OSV1	P_60_MSV	P_60_OSV1 	P_66_MSV 	P_66_OSV1	P_75_MSV	
		P_75_OSV1	P_80_MSV	P_80_OSV1 	P_95_MSV	P_95_OSV1	P_97_5_MSV	P_97_5_OSV1 	P_VAR_LIN 	
		P_VAR_S1	P_VAR_S2 	P_VAR_S3 	PB_KNOTS 	PB_REGRESSION	PB_SPECIF_VAL 	RCS 	REF_VAL_FA 	SD_MSV 	
		SD_OSV1	SP_VALUES_INDEX 	SP_VALUES1 	SP_VALUES2 	SP_VALUES3 	TEST_NAME 	TEXT 	TITRE_GPLOT 	
		TITRE_GPLOT_BF 	V_KNOT1 	V_KNOT1_MSV	V_KNOT1_OSV1	V_KNOT1_OSV2	V_KNOT1_OSV3	V_KNOT1_OSV4	
		V_KNOT1_OSV5	V_KNOT1_OSV6	V_KNOT1_OSV7	V_KNOT1_OSV8	V_KNOT1_OSV9	V_KNOT1_OSV10	V_KNOT2 	
		V_KNOT2_MSV 	V_KNOT2_OSV1 	V_KNOT2_OSV2	V_KNOT2_OSV3	V_KNOT2_OSV4	V_KNOT2_OSV5	V_KNOT2_OSV6	
		V_KNOT2_OSV7	V_KNOT2_OSV8	V_KNOT2_OSV9	V_KNOT2_OSV10	V_KNOT3 	V_KNOT3_MSV 	V_KNOT3_OSV1	
		V_KNOT3_OSV2	V_KNOT3_OSV3	V_KNOT3_OSV4	V_KNOT3_OSV5	V_KNOT3_OSV6	V_KNOT3_OSV7	V_KNOT3_OSV8	
		V_KNOT3_OSV9	V_KNOT3_OSV10	V_KNOT4 	V_KNOT4_MSV	V_KNOT4_OSV1	V_KNOT4_OSV2	V_KNOT4_OSV3	
		V_KNOT4_OSV4	V_KNOT4_OSV5	V_KNOT4_OSV6	V_KNOT4_OSV7	V_KNOT4_OSV8	V_KNOT4_OSV9	
		V_KNOT4_OSV10	V_KNOT5 	V_KNOT5_MSV	V_KNOT5_OSV1	V_KNOT5_OSV2	V_KNOT5_OSV3	V_KNOT5_OSV4	
		V_KNOT5_OSV5	V_KNOT5_OSV6	V_KNOT5_OSV7	V_KNOT5_OSV8	V_KNOT5_OSV9	V_KNOT5_OSV10	 	
		YAXIS_GPLOT MIN_Y MAX_Y  ;

%symdel 	AIC	FOR_RENAME 	I 	I_BIS 	INC	J 	K 	KNOT_INDEX 	KNOT1_MSV 	KNOT1_OSV1 	
			KNOT1_OSV2	KNOT1_OSV3	KNOT1_OSV4	KNOT1_OSV5	KNOT1_OSV6	KNOT1_OSV7	KNOT1_OSV8	KNOT1_OSV9	
			KNOT1_OSV10	KNOT2_MSV 	KNOT2_OSV1	KNOT2_OSV2	KNOT2_OSV3	KNOT2_OSV4	KNOT2_OSV5	KNOT2_OSV6	
			KNOT2_OSV7	KNOT2_OSV8	KNOT2_OSV9	KNOT2_OSV10	KNOT3_MSV 	KNOT3_OSV1 	KNOT3_OSV2	KNOT3_OSV3	KNOT3_OSV4	
			KNOT3_OSV5	KNOT3_OSV6	KNOT3_OSV7	KNOT3_OSV8	KNOT3_OSV9	KNOT3_OSV10	KNOT4_MSV 	KNOT4_OSV1	KNOT4_OSV2	
			KNOT4_OSV3	KNOT4_OSV4	KNOT4_OSV5	KNOT4_OSV6	KNOT4_OSV7	KNOT4_OSV8	KNOT4_OSV9	KNOT4_OSV10	KNOT5_MSV 	
			KNOT5_OSV1	KNOT5_OSV2	KNOT5_OSV3	KNOT5_OSV4	KNOT5_OSV5	KNOT5_OSV6	KNOT5_OSV7	KNOT5_OSV8	KNOT5_OSV9	
			KNOT5_OSV10	LEGEND 	LOG_L 	MAX_XAXIS_FA 	MEAN_MSV 	MEAN_OSV1 	MEAN_OSV2	MEAN_OSV3	
			MEAN_OSV4	MEAN_OSV5	MEAN_OSV6	MEAN_OSV7	MEAN_OSV8	MEAN_OSV9	MEAN_OSV10	MESSAGE_ADD1	
			MESSAGE_ADD2	MESSAGE_REF_VAL	MIN_XAXIS_FA 	MISTAKE 	NAME_MSV0 	NAME_MSV1 	NAME_MSV2 	NAME_MSV3 	
			NAME_OSV_CURR 	NAME10 	NAME11 	NB_CLASSES_BYFACT 	NB_KNOTS_CURR	NB_KNOTS_MINUS1 	NB_KNOTS_MSV 	
			NB_KNOTS_OSV1 	NB_KNOTS_OSV2	NB_KNOTS_OSV3	NB_KNOTS_OSV4	NB_KNOTS_OSV5	NB_KNOTS_OSV6	
			NB_KNOTS_OSV7	NB_KNOTS_OSV8	NB_KNOTS_OSV9	NB_KNOTS_OSV10	NB_OSV 	NB_PARAMETERS 	NB_SP_VALUES 	
			NO_CL_BYFACT 	OSV_CURR 	OUTFILE_FOR_ANALYSIS 	P_0_MSV	P_0_OSV1	P_100_MSV	P_100_OSV1	P_2_5_MSV 	
			P_2_5_OSV1 	P_20_MSV 	P_20_OSV1	P_25_MSV	P_25_OSV1	P_33_MSV 	P_33_OSV1	P_40_MSV	P_40_OSV1	
			P_5_MSV	P_5_OSV1 	P_50_MSV	P_50_OSV1	P_60_MSV	P_60_OSV1 	P_66_MSV 	P_66_OSV1	P_75_MSV	
			P_75_OSV1	P_80_MSV	P_80_OSV1 	P_95_MSV	P_95_OSV1	P_97_5_MSV	P_97_5_OSV1 	P_VAR_LIN 	
			P_VAR_S1	P_VAR_S2 	P_VAR_S3 	PB_KNOTS 	PB_REGRESSION	PB_SPECIF_VAL 	RCS 	REF_VAL_FA 	SD_MSV 	
			SD_OSV1	SP_VALUES_INDEX 	SP_VALUES1 	SP_VALUES2 	SP_VALUES3 	TEST_NAME 	TEXT 	TITRE_GPLOT 	
			TITRE_GPLOT_BF 	V_KNOT1 	V_KNOT1_MSV	V_KNOT1_OSV1	V_KNOT1_OSV2	V_KNOT1_OSV3	V_KNOT1_OSV4	
			V_KNOT1_OSV5	V_KNOT1_OSV6	V_KNOT1_OSV7	V_KNOT1_OSV8	V_KNOT1_OSV9	V_KNOT1_OSV10	V_KNOT2 	
			V_KNOT2_MSV 	V_KNOT2_OSV1 	V_KNOT2_OSV2	V_KNOT2_OSV3	V_KNOT2_OSV4	V_KNOT2_OSV5	V_KNOT2_OSV6	
			V_KNOT2_OSV7	V_KNOT2_OSV8	V_KNOT2_OSV9	V_KNOT2_OSV10	V_KNOT3 	V_KNOT3_MSV 	V_KNOT3_OSV1	
			V_KNOT3_OSV2	V_KNOT3_OSV3	V_KNOT3_OSV4	V_KNOT3_OSV5	V_KNOT3_OSV6	V_KNOT3_OSV7	V_KNOT3_OSV8	
			V_KNOT3_OSV9	V_KNOT3_OSV10	V_KNOT4 	V_KNOT4_MSV	V_KNOT4_OSV1	V_KNOT4_OSV2	V_KNOT4_OSV3	
			V_KNOT4_OSV4	V_KNOT4_OSV5	V_KNOT4_OSV6	V_KNOT4_OSV7	V_KNOT4_OSV8	V_KNOT4_OSV9	
			V_KNOT4_OSV10	V_KNOT5 	V_KNOT5_MSV	V_KNOT5_OSV1	V_KNOT5_OSV2	V_KNOT5_OSV3	V_KNOT5_OSV4	
			V_KNOT5_OSV5	V_KNOT5_OSV6	V_KNOT5_OSV7	V_KNOT5_OSV8	V_KNOT5_OSV9	V_KNOT5_OSV10		
			YAXIS_GPLOT MIN_Y MAX_Y;

%let mistake = 0;

/* ------------------------------------------- */
/* Number of other continuous spline variables */
/* ------------------------------------------- */

%let nb_osv = 0;
%do i = 1 %to 10;
%if &&oth_spline_var&i. ne %then %let nb_osv = %eval(&nb_osv. + 1);
%end;

/* ---------------- */
/* Knots definition */
/* ---------------- */


	/* For the main spline variable */

%let nb_knots_msv = %words(&knots_msv.);

data _null_;
nb_knots = &nb_knots_msv.;
if 0 < nb_knots < 3 then do;
call symput("mistake",1);
call symput("message"," The number of knots for &main_spline_var. is < 3 (it must be 3-5) ");
end;
if nb_knots > 5 then do;
call symput("mistake",1);
call symput("message"," The number of knots for &main_spline_var. is > 5 (it must be 3-5) ");
end;
run;

	%do knot_index = 1 %to &nb_knots_msv.;
	%let knot&knot_index._msv = %scan(&knots_msv.,&knot_index.,%STR( ));
	%end;

	/* For the other spline variables */

%do i = 1 %to &nb_osv.;

%let nb_knots_osv&i. = %words(&&knots_osv&i.);

data _null_;
nb_knots = &&nb_knots_osv&i.;
if 0 < nb_knots < 3 then do;
call symput("mistake",1);
call symput("message"," The number of knots for &&oth_spline_var&i. is < 3 (it must be 3-5) ");
end;
if nb_knots > 5 then do;
call symput("mistake",1);
call symput("message"," The number of knots for &&oth_spline_var&i. is > 5 (it must be 3-5) ");
end;
run;

	%do knot_index = 1 %to &&nb_knots_osv&i.;
	%let knot&knot_index._osv&i. = %scan(&&knots_osv&i.,&knot_index.,%STR( ));
	%put;
	%end;
%end;

%if %length(%str(&dir_data.)) > 0 %then %do;
libname _RCS "&dir_data.";
%let RCS = _RCS;
%end;

%if %length(%str(&dir_data.)) = 0 %then %do;
%let RCS = Work;
%end;

%if %length(%str(&exp_beta.)) = 0 %then %let exp_beta = 0;
%if %length(%str(&print_OR_HR.)) = 0 %then %let print_OR_HR = 0;
%if %length(%str(&histogram.)) = 0 %then %let histogram = 0;
%if %length(%str(&display_knots.)) = 0 %then %let display_knots = 1;
%if %length(%str(&X_ref_line.)) = 0 %then %let X_ref_line = 0;
%if %length(%str(&Y_ref_line.)) = 0 %then %let Y_ref_line = 0;
%if %length(%str(&no_graph.)) = 0 %then %let no_graph = 0;
%if %length(%str(&no_title.)) = 0 %then %let no_title = 0;
%if %length(%str(&no_label_X.)) = 0 %then %let no_label_X = 0;
%if %length(%str(&no_label_Y.)) = 0 %then %let no_label_Y = 0;
%if %length(%str(&no_legend.)) = 0 %then %let no_legend = 0;
%if %length(%str(&print_covar_mat.)) = 0 %then %let print_covar_mat = 0;
%if %length(%str(&missing_BF.)) = 0 %then %let missing_BF = 0;
%if %length(%str(&print_obs.)) = 0 %then %let print_obs = 0;
%if %length(%str(&width_band_obs.)) = 0 %then %let width_band_obs = 3;
%if %length(%str(&export_curves.)) = 0 %then %let export_curves = 0;

%if &AVK_msv. = %then %let AVK_msv = 0;
%do i = 1 %to 10;
%if &&AVK_osv&i. = %then %let AVK_osv&i. = 0;
%end;

%if %length(%str(&no_delete_files.)) = 0 %then %let no_delete_files = 0;
%if %length(%str(&subject_var.)) = 0 %then %let test_name = chisq;
%if %length(%str(&subject_var.)) > 0 %then %let test_name = Z;
%if %sysfunc(lowcase(%str(&typ_reg.))) = lin %then %do;%let for_rename = Differences;%end;
	%if %sysfunc(lowcase(%str(&typ_reg.))) = log %then %do;
	%if &exp_beta. = 0 %then %do;%let for_rename = Ln_ORs;%end;
	%if &exp_beta. = 1 %then %do;%let for_rename = ORs;%end;
	%end;
	%if %sysfunc(lowcase(%str(&typ_reg.))) = cox %then %do;
	%if &exp_beta. = 0 %then %do;%let for_rename = Ln_HRs;%end;
	%if &exp_beta. = 1 %then %do;%let for_rename = HRs;%end;
	%end;

%if &AVK_msv. = 0 %then %do;
	%if %length(%str(&knots_msv.)) = 0 %then %do;
	%let nb_knots_msv = 3;
	%let knot1_msv = 5;
	%let knot2_msv = 50;
	%let knot3_msv = 95;
	%end;
%end;

%do i = 1 %to &nb_osv.;

	%if &&AVK_osv&i. = 0 %then %do;
		%if %length(%str(&&knots_osv&i.)) = 0 %then %do;
		%let nb_knots_osv&i. = 3;
		%let knot1_osv&i. = 5;
		%let knot2_osv&i. = 50;
		%let knot3_osv&i. = 95;
		%end;
	%end;

%end;

/* ------------------------------- */
/* Creation of the infile datafile */
/* ------------------------------- */

data  &infile._II;
set &RCS..&infile.;
where &where.;
run ;

/* ------------------- */
/* Mistakes definition */
/* ------------------- */

%let message_add1 = ;
%let message_add2 = ;

%if %length(%str(&main_spline_var.)) > 12 %then %do;

data _null_;
old_name = "&main_spline_var.";
new_name = substr("&main_spline_var.",1,11 )||"_" ;
call symput ('old_name',old_name);
call symput ('new_name',new_name);
run ;

data  &infile._II;
set   &infile._II;
&new_name. = &main_spline_var.;
run ;

%let main_spline_var = &new_name.;
%let message_add1 = Note: the &old_name. variable contained more than 12 characters. It has therefore been truncated and its new name in the analyses is: &main_spline_var.;

%end;

%do i = 1 %to &nb_osv.;

	%if %length(%str(&&oth_spline_var&i.)) > 12 %then %do;

	data _null_;
	old_name = "&&oth_spline_var&i.";
	new_name = substr("&&oth_spline_var&i.",1,11 )||"_" ;
	call symput ('old_name',old_name);
	call symput ('new_name',new_name);
	run ;

	data  &infile._II;
	set   &infile._II;
	&new_name. = &&oth_spline_var&i.;
	run ;

	%let oth_spline_var&i. = &new_name.;
	%let message_add1 = Note: the &old_name. variable contained more than 12 characters. It has therefore been truncated and its new name in the analyses is: &&oth_spline_var&i.;

	%end;

%end;

%if %length(%str(&dep_var.)) > 16 %then %do;

data _null_;
set &infile._II;
old_name_dep_var = "&dep_var.";
new_name_dep_var = substr("&dep_var.",1,15 )||"_" ;
call symput ('old_name_dep_var',old_name_dep_var);
call symput ('new_name_dep_var',new_name_dep_var);
run ;

data  &infile._II;
set   &infile._II;
&new_name_dep_var. = &dep_var.;
run ;

%let dep_var = &new_name_dep_var.;
%let message_add2 = Note: the &old_name_dep_var. variable contained more than 16 characters. It has therefore been truncated and its new name in the analyses is: &dep_var.;

%end;

%if &AVK_msv. = 0 %then %do;
	%do j = 1 %to &nb_knots_msv.;
	data _null_;
	knot_curr = &&knot&j._msv;
		if knot_curr <= 0 then do;
		call symput("mistake",1);
		call symput("message","  The &j.th knot of &main_spline_var. has a negative value (a percentile must be > 0) ");
		end;
		if round(knot_curr,1) ne knot_curr then do;
		call symput("mistake",1);
		call symput("message","  A precentile must be an integer (see the value of the &j.th knot for &main_spline_var.) ");
		end;
	run;
	%end;
%end;

%do i = 1 %to &nb_osv.;

	%if &&AVK_osv&i. = 0 %then %do;
		%do j = 1 %to &&nb_knots_osv&i.;
		data _null_;
		knot_curr = &&knot&j._osv&i.;
			if knot_curr <= 0 then do;
			call symput("mistake",1);
			call symput("message","  The &j.th knot of &&oth_spline_var&i. has a negative value (a percentile must be > 0) ");
			end;
			if round(knot_curr,1) ne knot_curr then do;
			call symput("mistake",1);
			call symput("message","  A precentile must be an integer (see the value of the &j.th knot for &&oth_spline_var&i.) ");
			end;
		run;
		%end;
	%end;

%end;

%if %length(%str(&infile.)) = 0 %then %do;
%let mistake = 1;
%let message = " You forgot to specify the name of the original SAS datafile (INFILE parameter) ";
%end;

%if %length(%str(&subject_var.)) > 0 and %length(%str(&work_corr_matrix.)) = 0 %then %do;
%let mistake = 1;
%let message = " You forgot to specify the working correlation matrix structure for your GEE model (WORK_CORR_MATRIX parameter) ";
%end;

%if %length(%str(&main_spline_var.)) = 0 %then %do;
%let mistake = 1;
%let message = " You forgot to specify the name of the main continuous exposure (MAIN_SPLINE_VAR parameter) ";
%end;

%if &AVK_msv. = 1 %then %do;
	%if	&nb_knots_msv. = 0 %then %do;;
	%let mistake = 1;
	%let message = " By specifying AVK = 1 for &main_spline_var., you must provide the values of at least 3 knots ";
	%end;
%end;

%if &export_curves. = 1 %then %do;
	%if	%length(%str(&dir_export_curves.)) = 0 %then %do;;
	%let mistake = 1;
	%let message = " You forgot to specify the directory where you want to export the txt file for curves ";
	%end;
%end;

data _null_;
if &AVK_msv. not in (0,1) then do;
call symput("mistake",1);
call symput("message"," The AVK parameter for &main_spline_var. must value 0 or 1 (you typed: '&AVK_msv.')");
end;
run;

%do i = 1 %to &nb_osv.;

	%if &&AVK_osv&i. = 1 %then %do;
		%if	&&nb_knots_osv&i. = 0 %then %do;;
		%let mistake = 1;
		%let message = " By specifying AVK = 1 for &&oth_spline_var&i., you must provide the values of at least 3 knots ";
		%end;
	%end;

	data _null_;
	if &&AVK_osv&i. not in (0,1) then do;
	call symput("mistake",1);
	call symput("message"," The AVK parameter for &&oth_spline_var&i. must value 0 or 1 (you typed: '&&AVK_osv&i.')");
	end;
	run;

%end;

%if &typ_reg. = lin and &exp_beta. = 1 %then %do;
%let exp_beta = 0;
%end;

%if %length(%str(&typ_reg.)) > 0 %then %do;
data _null_;
typ_reg = "%sysfunc(lowcase(%str(&typ_reg.)))";
call symput ("typ_reg",typ_reg);
if typ_reg not in ("lin","log","cox") then do;
call symput("mistake",1);
call symput("message"," You made a mistake when you specified the TYP_REG parameter (you typed: '&typ_reg.') ");
end;
run;
%end;

%if %length(%str(&typ_reg.)) > 0 %then %do;

	%if %length(%str(&dep_var.)) = 0 %then %do;
	%let mistake = 1;
	%let message = " You forgot to specify the name of the dependent variable for the regression (DEP_VAR parameter) ";
	%end;

	%if &typ_reg. = cox %then %do;

		%if %length(%str(&surv_time_var.)) = 0 %then %do;
		%let mistake = 1;
		%let message = " You forgot to specify the name of the survival time variable (SURV_TIME_VAR parameter) ";
		%end;

	%end;

%end;
%if &mistake. = 1 %then %do;
%put ;
%put ;
%put *		+------------------+;
%put *		| **** ERROR! **** |;
%put *		+------------------+;
%put ;
%put &message.;
%end;

/* -------------------------- */
/* End of mistakes definition */
/* -------------------------- */

	/* Definition of classes of by_factor variable */

%if %length(%str(&by_factor.)) > 0 %then %do;
Proc freq data = &infile._II noprint;
tables &by_factor. / out= data_by_factor;
run;

data data_by_factor;
set data_by_factor;
	%if &missing_BF. = 0 %then %do;
	if &by_factor. = . then delete;
	%end;
	%if &missing_BF. = 1 %then %do;
	if &by_factor. = . then &by_factor. = -99;
	%end;
run;

proc sort data = data_by_factor;
by &by_factor.;
run;

data _null_;
set data_by_factor;
by &by_factor.;
if last.&by_factor. then call symput("nb_classes_byfact",_n_);
run;

%let nb_classes_byfact = %trim(%left(&nb_classes_byfact.));

data _null_;
set data_by_factor;
	%do i = 1 %to &nb_classes_byfact.;
	if _n_ = &i. then call symput ("VBF_&i.",%trim(%left(&by_factor.)));
	%end;
run;

	%do i = 1 %to &nb_classes_byfact.;
	%let val_by_factor_&i. = %trim(%left(&&VBF_&i.));
	%end;

data &infile._II;
set &infile._II;
if &by_factor. = . then &by_factor. = -99;
run;

%end;


%if %length(%str(&by_factor.)) = 0 %then %let nb_classes_byfact = 1;

%do no_cl_byfact = 1 %to &nb_classes_byfact.;	/* Beginning of BY_FACTOR loop */

	%if %length(%str(&by_factor.)) = 0 %then %do;
	data &infile._for_analysis;
	set &infile._II;
	run;

		%if %length(%str(&outfile.)) = 0 %then %do;
		%let outfile_for_analysis = &infile._RCS;
		%end;

		%if %length(%str(&outfile.)) > 0 %then %do;
		%let outfile_for_analysis = &outfile.;
		%end;

	%end;

	%if %length(%str(&by_factor.)) > 0 %then %do;
	data &infile._for_analysis;
	set &infile._II;
	where &by_factor. = &&val_by_factor_&no_cl_byfact.;
	run;

		%if %length(%str(&outfile.)) = 0 %then %do;
		%let outfile_for_analysis = &infile._RCS_&no_cl_byfact.;
		%end;

		%if %length(%str(&outfile.)) > 0 %then %do;
		%let outfile_for_analysis = &outfile._&no_cl_byfact.;
		%end;

	%end;

%if &mistake. = 0 %then %do;	/* Beginning of mistake-free loop #1 */

%if &AVK_msv. = 0 %then %do;

	%if &nb_knots_msv. = 3 %then %do;

	proc univariate data = &infile._for_analysis noprint;
	var &main_spline_var.;
	output out = out_univariate	n = nbre pctlpre=P_ 
						pctlpts=&knot1_msv.,&knot2_msv.,&knot3_msv.;
	run;

	data Defknots_from_univ;
	set out_univariate;
	knot1 = P_&knot1_msv.;
	knot2 = P_&knot2_msv.;
	knot3 = P_&knot3_msv.;
	knot4 = .;
	knot5 = .;
	run;

	%end;

	%if &nb_knots_msv. = 4 %then %do;

	proc univariate data = &infile._for_analysis noprint;
	var &main_spline_var.;
	output out = out_univariate	n = nbre pctlpre=P_ 
						pctlpts=&knot1_msv.,&knot2_msv.,&knot3_msv.,&knot4_msv.;
	run;

	data Defknots_from_univ;
	set out_univariate;
	knot1 = P_&knot1_msv.;
	knot2 = P_&knot2_msv.;
	knot3 = P_&knot3_msv.;
	knot4 = P_&knot4_msv.;
	knot5 = .;
	run;

	%end;

	%if &nb_knots_msv. = 5 %then %do;

	proc univariate data = &infile._for_analysis noprint;
	var &main_spline_var.;
	output out = out_univariate	n = nbre pctlpre=P_ 
						pctlpts=&knot1_msv.,&knot2_msv.,&knot3_msv.,&knot4_msv.,&knot5_msv.;
	run;

	data Defknots_from_univ;
	set out_univariate;
	knot1 = P_&knot1_msv.;
	knot2 = P_&knot2_msv.;
	knot3 = P_&knot3_msv.;
	knot4 = P_&knot4_msv.;
	knot5 = P_&knot5_msv.;
	run;

	%end;

data _null_;
set Defknots_from_univ;
call symput ("V_knot1_msv",knot1);
call symput ("V_knot2_msv",knot2);
call symput ("V_knot3_msv",knot3);
call symput ("V_knot4_msv",knot4);
call symput ("V_knot5_msv",knot5);
run;

%end;

%if &AVK_msv. = 1 %then %do;

	%do j = 1 %to &nb_knots_msv.;
	data _null_;
	knot&j._msv = &&knot&j._msv;
	call symput ("V_knot&j._msv",knot&j._msv);
	run;
	%end;

%end;

%do i = 1 %to &nb_osv.;

	%if &&AVK_osv&i. = 0 %then %do;

		%if &&nb_knots_osv&i. = 3 %then %do;

		proc univariate data = &infile._for_analysis noprint;
		var &&oth_spline_var&i.;
		output out = out_univariate	n = nbre pctlpre=P_ 
						pctlpts=&&knot1_osv&i.,&&knot2_osv&i.,&&knot3_osv&i.;
		run;

		data Defknots_from_univ;
		set out_univariate;
		knot1 = P_&&knot1_osv&i.;
		knot2 = P_&&knot2_osv&i.;
		knot3 = P_&&knot3_osv&i.;
		knot4 = .;
		knot5 = .;
		run;

		%end;

		%if &&nb_knots_osv&i. = 4 %then %do;

		proc univariate data = &infile._for_analysis noprint;
		var &&oth_spline_var&i.;
		output out = out_univariate	n = nbre pctlpre=P_ 
						pctlpts=&&knot1_osv&i.,&&knot2_osv&i.,&&knot3_osv&i.,&&knot4_osv&i.;
		run;

		data Defknots_from_univ;
		set out_univariate;
		knot1 = P_&&knot1_osv&i.;
		knot2 = P_&&knot2_osv&i.;
		knot3 = P_&&knot3_osv&i.;
		knot4 = P_&&knot4_osv&i.;
		knot5 = .;
		run;

		%end;

		%if &&nb_knots_osv&i. = 5 %then %do;

		proc univariate data = &infile._for_analysis noprint;
		var &&oth_spline_var&i.;
		output out = out_univariate	n = nbre pctlpre=P_ 
						pctlpts=&&knot1_osv&i.,&&knot2_osv&i.,&&knot3_osv&i.,&&knot4_osv&i.,&&knot5_osv&i.;
		run;

		data Defknots_from_univ;
		set out_univariate;
		knot1 = P_&&knot1_osv&i.;
		knot2 = P_&&knot2_osv&i.;
		knot3 = P_&&knot3_osv&i.;
		knot4 = P_&&knot4_osv&i.;
		knot5 = P_&&knot5_osv&i.;
		run;

		%end;

	data _null_;
	set Defknots_from_univ;
	call symput ("V_knot1_osv&i.",knot1);
	call symput ("V_knot2_osv&i.",knot2);
	call symput ("V_knot3_osv&i.",knot3);
	call symput ("V_knot4_osv&i.",knot4);
	call symput ("V_knot5_osv&i.",knot5);
	run;

	%end;

	%if &&AVK_osv&i. = 1 %then %do;

		%do j = 1 %to &&nb_knots_osv&i.;
		data _null_;
		knot&j._osv&i. = &&knot&j._osv&i.;
		call symput ("V_knot&j._osv&i.",knot&j._osv&i.);
		run;
		%end;

	%end;

%end;

%put;
%put;
%put;
%put;
%if %length(%str(&by_factor.)) = 0 %then %put ---------------------------- Miscellaneous information ----------------------------;
%if %length(%str(&by_factor.)) > 0 %then %put ------------------------ Miscellaneous information when &by_factor. = &&val_by_factor_&no_cl_byfact. ------------------------;
%put;
%put;
%if &AVK_msv. = 0 %then %do;
%put Values of the &nb_knots_msv. knots of &main_spline_var. (choice according to the percentiles):;
	%do j = 1 %to &nb_knots_msv.;
	%put Knot #&j.: %trim(%left(&&V_knot&j._msv)) 		(&&knot&j._msv th percentile);
	%end;
%end;
%if &AVK_msv. = 1 %then %do;
%put Values of the &nb_knots_msv. knots of &main_spline_var. (chosen by the user):;
	%do j = 1 %to &nb_knots_msv.;
	%put Knot #&j.: %trim(%left(&&V_knot&j._msv));
	%end;
%end;

%do i = 1 %to &nb_osv.;
%put;
%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
%put;
	%if &&AVK_osv&i. = 0 %then %do;
	%put Values of the &&nb_knots_osv&i. knots of &&oth_spline_var&i. (choice according to the percentiles):;
		%do j = 1 %to &&nb_knots_osv&i.;
		%put Knot #&j.: %trim(%left(&&V_knot&j._osv&i.)) 		(&&knot&j._osv&i. th percentile);
		%end;
	%end;
	%if &&AVK_osv&i. = 1 %then %do;
	%put Values of the &&nb_knots_osv&i. knots of &&oth_spline_var&i. (chosen by the user):;
		%do j = 1 %to &&nb_knots_osv&i.;
		%put Knot #&j.: %trim(%left(&&V_knot&j._osv&i.));
		%end;
	%end;
%end;

proc univariate data = &infile._for_analysis noprint;
var &main_spline_var.;
output out = out_distrib   n = nbre mean = mean std = sd
 	                          pctlpre=P_ pctlpts=0,2.5,5,20,25,33,40,50,60,66,75,80,95,97.5,100;
run;

data _null_;
set out_distrib;
call symput ("P_0_msv",P_0);
call symput ("P_2_5_msv",P_2_5);
call symput ("P_5_msv",P_5);
call symput ("P_20_msv",P_20);
call symput ("P_25_msv",P_25);
call symput ("P_33_msv",P_33);
call symput ("P_40_msv",P_40);
call symput ("P_50_msv",P_50);
call symput ("P_60_msv",P_60);
call symput ("P_66_msv",P_66);
call symput ("P_75_msv",P_75);
call symput ("P_80_msv",P_80);
call symput ("P_95_msv",P_95);
call symput ("P_97_5_msv",P_97_5);
call symput ("P_100_msv",P_100);
call symput ("Mean_msv",mean);
call symput ("SD_msv",sd);
run;

%let pb_knots = 0;
%let pb_specif_val = 0;

%do i = 1 %to %eval(&nb_knots_msv. - 1);
%let j = %eval(&i. + 1);
data _null_;
if &&V_knot&i._msv = &&V_knot&j._msv then call symput("pb_knots",1);
if &&V_knot&i._msv > &&V_knot&j._msv then call symput("pb_knots",3);
run;
%end;

%if &AVK_msv. = 1 %then %do;
	%do i = 1 %to &nb_knots_msv.;
	%let j = %eval(&i. + 1);
	data _null_;
	if &&V_knot&i._msv > &P_100_msv. or &&V_knot&i._msv < &P_0_msv. then call symput("pb_knots",2);
	run;
	%end;
%end;

%if &pb_knots. = 1 %then %do;
%let mistake = 1;
%let message = " At least 2 adjacents knots are equal for &main_spline_var. -> please check the location of your knots for this variable";
%end;

%if &pb_knots. = 2 %then %do;
%let mistake = 1;
%let message = " One of the &nb_knots_msv. knots for &main_spline_var. is greater than the maximum of &main_spline_var. or lower of the minimum of &main_spline_var. ";
%end;

%if &pb_knots. = 3 %then %do;
%let mistake = 1;
%let message = " Knots for &main_spline_var. have not been listed in ascendant order ";
%end;

%if %length(%str(&specif_val.)) > 0 and %length(%str(&typ_reg.)) > 0 %then %do;

%let nb_sp_values = %words(&specif_val.);
	%do sp_values_index = 1 %to &nb_sp_values.;
	%let sp_values&sp_values_index. = %scan(&specif_val.,&sp_values_index.,%STR( ));
	data _null_;
	if (&&sp_values&sp_values_index. > &P_100_msv. or &&sp_values&sp_values_index. < &P_0_msv.) then call symput("pb_specif_val",1);
	run;
	%end;
%end;

%if &pb_specif_val. = 1 %then %do;
%let mistake = 1;
%let message = " One of the value(s) that you assigned for the SPECIF_VAL parameter for &main_spline_var. is greater than the maximum of &main_spline_var. or lower of the minimum of &main_spline_var.";
%end;

%do k = 1 %to &nb_osv.;

	proc univariate data = &infile._for_analysis noprint;
	var &&oth_spline_var&k.;
	output out = out_distrib   n = nbre mean = mean std = sd
 	                          pctlpre=P_ pctlpts=0,2.5,5,20,25,33,40,50,60,66,75,80,95,97.5,100;
	run;

	data _null_;
	set out_distrib;
	call symput ("P_0_osv&k.",P_0);
	call symput ("P_2_5_osv&k.",P_2_5);
	call symput ("P_5_osv&k.",P_5);
	call symput ("P_20_osv&k.",P_20);
	call symput ("P_25_osv&k.",P_25);
	call symput ("P_33_osv&k.",P_33);
	call symput ("P_40_osv&k.",P_40);
	call symput ("P_50_osv&k.",P_50);
	call symput ("P_60_osv&k.",P_60);
	call symput ("P_66_osv&k.",P_66);
	call symput ("P_75_osv&k.",P_75);
	call symput ("P_80_osv&k.",P_80);
	call symput ("P_95_osv&k.",P_95);
	call symput ("P_97_5_osv&k.",P_97_5);
	call symput ("P_100_osv&k.",P_100);
	call symput ("Mean_osv&k.",mean);
	call symput ("SD_osv&k.",sd);
	run;

	%let pb_knots = 0;
	%let pb_specif_val = 0;

	%do i = 1 %to %eval(&&nb_knots_osv&k. - 1);
	%let j = %eval(&i. + 1);
	data _null_;
	if &&V_knot&i._osv&k. = &&V_knot&j._osv&k. then call symput("pb_knots",1);
	if &&V_knot&i._osv&k. > &&V_knot&j._osv&k. then call symput("pb_knots",3);
	run;
	%end;

	%if &&AVK_osv&k. = 1 %then %do;
		%do i = 1 %to &&nb_knots_osv&k.;
		%let j = %eval(&i. + 1);
		data _null_;
		if &&V_knot&i._osv&k. > &&P_100_osv&k. or &&V_knot&i._osv&k. < &&P_0_osv&k. then call symput("pb_knots",2);
		run;
		%end;
	%end;

	%if &pb_knots. = 1 %then %do;
	%let mistake = 1;
	%let message = " At least 2 adjacents knots are equal for &&oth_spline_var&k. -> please check the location of your knots for this variable";
	%end;

	%if &pb_knots. = 2 %then %do;
	%let mistake = 1;
	%let message = "- One of the &&nb_knots_osv&k. knots for &&oth_spline_var&k. is greater than the maximum of &&oth_spline_var&k. or lower of the minimum of &&oth_spline_var&k. ";
	%end;

	%if &pb_knots. = 3 %then %do;
	%let mistake = 1;
	%let message = " Knots for &&oth_spline_var&k. have not been listed in ascendant order ";
	%end;

%end;

%if &mistake. = 1 %then %do;

	%if %length(%str(&by_factor.)) > 0 %then %do;
	%put ;
	%put ;
	%put ;
	%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
	%put For &by_factor. = &&val_by_factor_&no_cl_byfact.;
	%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
	%end;

%put ;
%put ;
%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
%put ;
%put Distribution of &main_spline_var. to check errors:;
%put ;
%put Minimum           =  %trim(%left(&&P_0_msv.));
%put 2.5th percentile  =  %trim(%left(&&P_2_5_msv.));
%put 5th   percentile  =  %trim(%left(&&P_5_msv.));
%put 25th percentile   =  %trim(%left(&&P_25_msv.));
%put Median            =  %trim(%left(&&P_50_msv.));
%put 75th percentile   =  %trim(%left(&&P_75_msv.));
%put 95th percentile   =  %trim(%left(&&P_95_msv.));
%put 97.5th percentile =  %trim(%left(&&P_97_5_msv.));
%put Maximum           =  %trim(%left(&&P_100_msv.));
%put ;

	%do i = 1 %to &nb_osv.;
	%put ;
	%put ;
	%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
	%put ;
	%put Distribution of &&oth_spline_var&i. to check errors:;
	%put ;
	%put Minimum           =  %trim(%left(&&P_0_osv&i.));
	%put 2.5th percentile  =  %trim(%left(&&P_2_5_osv&i.));
	%put 5th   percentile  =  %trim(%left(&&P_5_osv&i.));
	%put 25th percentile   =  %trim(%left(&&P_25_osv&i.));
	%put Median            =  %trim(%left(&&P_50_osv&i.));
	%put 75th percentile   =  %trim(%left(&&P_75_osv&i.));
	%put 95th percentile   =  %trim(%left(&&P_95_osv&i.));
	%put 97.5th percentile =  %trim(%left(&&P_97_5_osv&i.));
	%put Maximum           =  %trim(%left(&&P_100_osv&i.));
	%put ;
	%end;

%put ;
%put ;
%put &message.;
%put ;
%put *		+--------------------------------------+;
%put *		| **** ERROR (See message above)! **** |;
%put *		+--------------------------------------+;
%put ;
%end;

%end;	/* End of mistake-free loop #1 */



%if &mistake. = 0 %then %do;		/* Beginning of mistake-free loop #2 */

%if %length(%str(&min_Xaxis.)) > 0 %then %let min_Xaxis_FA = &min_Xaxis.;
%if %length(%str(&min_Xaxis.)) = 0 %then %do;
%let min_Xaxis_FA = &P_0_msv.;
%end;

%if %length(%str(&max_Xaxis.)) > 0 %then %let max_Xaxis_FA = &max_Xaxis.;
%if %length(%str(&max_Xaxis.)) = 0 %then %do;
%let max_Xaxis_FA = &P_100_msv.;
%end;

%let message_ref_val =;
%if %length(%str(&ref_val.)) > 0 and %length(%str(&typ_reg.)) > 0 %then %let ref_val_FA = &ref_val.;
%if %length(%str(&ref_val.)) = 0 and %length(%str(&typ_reg.)) > 0 %then %do;
%let ref_val_FA = &P_50_msv.;
%let message_ref_val = The reference value for &main_spline_var. has been assigned to its median (i.e., %trim(%left(&p_50_msv.)));
%end;

%if &histogram. = 1 %then %do;
%if %length(%str(&by_factor.)) = 0 %then %do;Title "Distribution of &main_spline_var.";%end;
%if %length(%str(&by_factor.)) > 0 %then %do;Title "Distribution of &main_spline_var. for &by_factor. = &&val_by_factor_&no_cl_byfact.";%end;
Proc univariate data  = &infile._for_analysis noprint ;
Histogram &main_spline_var. / CFILL = yellow ;
RUN  ;
%end;

title;

%do i = 0 %to %eval(&nb_knots_msv. - 2);
	%if &i. = 0 %then %do;
	data _null_;
	name = cat("&main_spline_var."!!"_RCS_lin");
	call symput("name_msv&i.",name);
	run;
	%end;
	%if &i. > 0 %then %do;
	data _null_;
	name = cat("&main_spline_var.","_RCS_S&i");
	call symput("name_msv&i.",name);
	run;
	%end;
%end;

%do j = 1 %to &nb_osv.;
	%do i = 0 %to %eval(&&nb_knots_osv&j. - 2);
		%if &i. = 0 %then %do;
		data _null_;
		name = cat("&&oth_spline_var&j."!!"_RCS_lin");
		call symput("name&j.&i.",name);
		run;
		%end;
		%if &i. > 0 %then %do;
		data _null_;
		name = cat("&&oth_spline_var&j.","_RCS_S&i");
		call symput("name&j.&i.",name);
		run;
		%end;
	%end;
%end;

/* -------------------------- */
/* Data creation with splines */
/* -------------------------- */

data &outfile_for_analysis.;
set &infile._for_analysis;

name_spline_msv0 = &main_spline_var.;

%let nb_knots_curr = &nb_knots_msv.;
%let nb_knots_minus1 = %eval(&nb_knots_msv. - 1);

%do i = 1 %to %eval(&nb_knots_msv. - 2);

name_spline_msv&i. =	(&main_spline_var. > &&V_knot&i._msv)*(&main_spline_var. - &&V_knot&i._msv)**3  -
						(&main_spline_var. > &&V_knot&nb_knots_minus1._msv)  * ((&&V_knot&nb_knots_curr._msv - &&V_knot&i._msv)/(&&V_knot&nb_knots_curr._msv - &&V_knot&nb_knots_minus1._msv)) * (&main_spline_var. - &&V_knot&nb_knots_minus1._msv)**3 +
						(&main_spline_var. > &&V_knot&nb_knots_curr._msv) * ((&&V_knot&nb_knots_minus1._msv - &&V_knot&i._msv)/(&&V_knot&nb_knots_curr._msv - &&V_knot&nb_knots_minus1._msv)) * (&main_spline_var. - &&V_knot&nb_knots_curr._msv)**3
						;
%end;

%do j = 1 %to &nb_osv.;

name_spline_&j.0 = &&oth_spline_var&j.;

%let nb_knots_curr = &&nb_knots_osv&j.;
%let nb_knots_minus1 = %eval(&&nb_knots_osv&j. - 1);

		%do i = 1 %to %eval(&&nb_knots_osv&j. - 2);

		name_spline_&j.&i. =	(&&oth_spline_var&j. > &&V_knot&i._osv&j.)*(&&oth_spline_var&j. - &&V_knot&i._osv&j.)**3  -
								(&&oth_spline_var&j. > &&V_knot&nb_knots_minus1._osv&j.)  * ((&&V_knot&nb_knots_curr._osv&j. - &&V_knot&i._osv&j.)/(&&V_knot&nb_knots_curr._osv&j. - &&V_knot&nb_knots_minus1._osv&j.)) * (&&oth_spline_var&j. - &&V_knot&nb_knots_minus1._osv&j.)**3 +
								(&&oth_spline_var&j. > &&V_knot&nb_knots_curr._osv&j.) * ((&&V_knot&nb_knots_minus1._osv&j. - &&V_knot&i._osv&j.)/(&&V_knot&nb_knots_curr._osv&j. - &&V_knot&nb_knots_minus1._osv&j.)) * (&&oth_spline_var&j. - &&V_knot&nb_knots_curr._osv&j.)**3
								;
		%end;

%end;


%put;
%do i = 0 %to %eval(&nb_knots_msv. - 2);
rename name_spline_msv&i. = &&name_msv&i.;
%end;

%do j = 1 %to &nb_osv.;
%put;
	%do i = 0 %to %eval(&&nb_knots_osv&j. - 2);
	rename name_spline_&j.&i. = &&name&j.&i.;
	%end;
%end;

run;

/* ------------------------------------- */
/* PUT instructions in the OUTPUT window */
/* ------------------------------------- */

data _null_;
file print ;
title1 " ";
%if %length(%str(&by_factor.)) = 0 %then %do;title2 "Distribution of the spline variable(s)";%end;
%if %length(%str(&by_factor.)) > 0 %then %do;title2 "Distribution of the spline variable(s) for &by_factor. = &&val_by_factor_&no_cl_byfact.";%end;
%if %length(%str(&by_factor.)) = 0 %then %do;title3 "--------------------------------------";%end;
%if %length(%str(&by_factor.)) > 0 %then %do;title3 "--------------------------------------------------------";%end;
put ;
put "------------------------------";
put ;
put "Distribution of &main_spline_var.:";
put ;
put "Minimum           =  %trim(%left(&P_0_msv.))";
put "2.5th percentile  =  %trim(%left(&P_2_5_msv.))";
put "5th   percentile  =  %trim(%left(&P_5_msv.))";
put "20th percentile   =  %trim(%left(&P_20_msv.))";
put "25th percentile   =  %trim(%left(&P_25_msv.))";
put "33th percentile   =  %trim(%left(&P_33_msv.))";
put "40th percentile   =  %trim(%left(&P_40_msv.))";
put "Median            =  %trim(%left(&P_50_msv.))";
put "Mean (SD)         =  %trim(%left(&mean_msv.)) (&SD_msv.)";
put "60th percentile   =  %trim(%left(&P_60_msv.))";
put "66th percentile   =  %trim(%left(&P_66_msv.))";
put "75th percentile   =  %trim(%left(&P_75_msv.))";
put "80th percentile   =  %trim(%left(&P_80_msv.))";
put "95th percentile   =  %trim(%left(&P_95_msv.))";
put "97.5th percentile =  %trim(%left(&P_97_5_msv.))";
put "Maximum           =  %trim(%left(&P_100_msv.))";
put ;
put "------------------------------";
put ;
%do i = 1 %to &nb_osv.;
put "Distribution of &&oth_spline_var&i.:";
put ;
put "Minimum           =  %trim(%left(&&P_0_osv&i.))";
put "2.5th percentile  =  %trim(%left(&&P_2_5_osv&i.))";
put "5th   percentile  =  %trim(%left(&&P_5_osv&i.))";
put "20th percentile   =  %trim(%left(&&P_20_osv&i.))";
put "25th percentile   =  %trim(%left(&&P_25_osv&i.))";
put "33th percentile   =  %trim(%left(&&P_33_osv&i.))";
put "40th percentile   =  %trim(%left(&&P_40_osv&i.))";
put "Median            =  %trim(%left(&&P_50_osv&i.))";
put "Mean (SD)         =  %trim(%left(&&mean_osv&i.)) (&&SD_osv&i.)";
put "60th percentile   =  %trim(%left(&&P_60_osv&i.))";
put "66th percentile   =  %trim(%left(&&P_66_osv&i.))";
put "75th percentile   =  %trim(%left(&&P_75_osv&i.))";
put "80th percentile   =  %trim(%left(&&P_80_osv&i.))";
put "95th percentile   =  %trim(%left(&&P_95_osv&i.))";
put "97.5th percentile =  %trim(%left(&&P_97_5_osv&i.))";
put "Maximum           =  %trim(%left(&&P_100_osv&i.))";
put ;
put "------------------------------";
%end;
put ;
put ;
run;

/* ----------- */
/* Regressions */
/* ----------- */

%let Pb_regression = 0;
%if %length(%str(&typ_reg.)) > 0 %then %do;

title ;

%let Titre_gplot = "Association between &dep_var. and &main_spline_var. using RCS with &nb_knots_msv. knots";

	%if %length(%str(&by_factor.)) = 0 %then %do;
	%let Titre_gplot_BF = " ";
	%end;
	%if %length(%str(&by_factor.)) > 0 %then %do;
	%let Titre_gplot_BF = "&by_factor. = &&val_by_factor_&no_cl_byfact.";
	%end;

	 /* ------------------- */
	 /* Logistic regression */
	 /* ------------------- */

	%if &typ_reg. = log %then %do;

	  	%if &exp_beta. = 0 %then %do;
	 	%let text = "Ln(OR) where the reference for &main_spline_var. is &ref_val_FA.";
	 	%let Yaxis_gplot = "Ln(OR) where the ref value for &main_spline_var. is &ref_val_FA.";
		%end;
		%if &exp_beta. = 1 %then %do;
		%let text = "OR where the reference for &main_spline_var. is &ref_val_FA.";
	 	%let Yaxis_gplot = "OR where the ref value for &main_spline_var. is &ref_val_FA.";
		%end;

	title1 " ";
	%if %length(%str(&by_factor.)) = 0 %then %do;title2 "   Logistic regression model including '&main_spline_var.' using a RCS function with &nb_knots_msv. knots";%end;
	%if %length(%str(&by_factor.)) > 0 %then %do;title2 "   Logistic regression model including '&main_spline_var.' using a RCS function with &nb_knots_msv. knots for &by_factor. = &&val_by_factor_&no_cl_byfact.";%end;
	title3 " ";
	title4 "             * - * - * - * - * - * - * - * - * - * - * - * - *";
	title5 " ";
 
	 	%if %length(%str(&subject_var.)) = 0 %then %do;
		ods output Modelfit = data_fit;
		ods output ParameterEstimates = data_parameters;
		ods output CovB = data_mat_cov;
		%if &print_covar_mat. = 0 %then %do;ods exclude CovB;%end;
		%end;

	 	%if %length(%str(&subject_var.)) > 0 %then %do;
		ods output GEEEmpPEst = data_parameters;
		ods output GEERCov =    data_mat_cov;
		%if &print_covar_mat. = 0 %then %do;ods exclude GEENCov CovB GEERCov;%end;
		%end;

	proc genmod data = &outfile_for_analysis. desc ;
	%if %length(%str(&subject_var.)) > 0 %then %do;class &subject_var.;%end;
	model &dep_var. = 	&main_spline_var._RCS_lin &main_spline_var._RCS_S1-&main_spline_var._RCS_S%eval(&nb_knots_msv. - 2)
						%do i = 1 %to &nb_osv.;
						%let osv_curr = &&oth_spline_var&i.;
						&osv_curr._RCS_lin &osv_curr._RCS_S1-&osv_curr._RCS_S%eval(&&nb_knots_osv&i. - 2) 
						%end;
						&adjust_var. / dist = binomial link = logit covb;
	
		output out = &outfile_for_analysis. predicted = Pred_&dep_var. lower = Lo_Pred_&dep_var. upper = Hi_Pred_&dep_var.;

	 	%if %length(%str(&subject_var.)) > 0 %then %do;
		repeated subject = &subject_var. / type = &work_corr_matrix. CovB;
		%end;

		%if &nb_knots_msv. = 3 %then %do;
		contrast "Overall_association" &main_spline_var._RCS_lin 1,&main_spline_var._RCS_S1 1 / wald;
		contrast "Non_lin_association" &main_spline_var._RCS_S1 1 / wald;
		%end;

		%if &nb_knots_msv. = 4 %then %do;
		contrast "Overall_association" &main_spline_var._RCS_lin 1,&main_spline_var._RCS_S1 1,&main_spline_var._RCS_S2 1 / wald;
		contrast "Non_lin_association" &main_spline_var._RCS_S1 1, &main_spline_var._RCS_S2 1 / wald;
		%end;

		%if &nb_knots_msv. = 5 %then %do;
		contrast "Overall_association" &main_spline_var._RCS_lin 1,&main_spline_var._RCS_S1 1,&main_spline_var._RCS_S2 1,&main_spline_var._RCS_S3 1 / wald;
		contrast "Non_lin_association" &main_spline_var._RCS_S1 1,&main_spline_var._RCS_S2 1,&main_spline_var._RCS_S3 1 / wald;
		%end;

	run;

	quit;
	ods output close;

 	%end;

	/* ----------*/
	/* Cox model */
	/* ----------*/

	%if &typ_reg. = cox %then %do;

	 	%if &exp_beta. = 0 %then %do;
	 	%let text = "Ln(HR) where the reference for &main_spline_var. is &ref_val_FA.";
	 	%let Yaxis_gplot = "Ln(HR) where the ref value for &main_spline_var. is &ref_val_FA.";
		%end;
		%if &exp_beta. = 1 %then %do;
		%let text = "HR where the reference for &main_spline_var. is &ref_val_FA.";
	 	%let Yaxis_gplot = "HR where the ref value for &main_spline_var. is &ref_val_FA.";
		%end;

	ods output ParameterEstimates = data_parameters;
	ods output FitStatistics = data_fit;
	ods output CovB = data_mat_cov;
	%if &print_covar_mat. = 0 %then %do;ods exclude CovB;%end;


	title1 " ";
	%if %length(%str(&by_factor.)) = 0 %then %do;title2 "   Cox model including '&main_spline_var.' using a RCS function with &nb_knots_msv. knots";%end;
	%if %length(%str(&by_factor.)) > 0 %then %do;title2 "   Cox model including '&main_spline_var.' using a RCS function with &nb_knots_msv. knots for &by_factor. = &&val_by_factor_&no_cl_byfact.";%end;
	title2 "   Cox model including '&main_spline_var.' using a RCS function with &nb_knots_msv. knots";
	title3 " ";
	title4 "             * - * - * - * - * - * - * - * - * - * - *";
	title5 " ";
	proc phreg data = &outfile_for_analysis.;
	model &surv_time_var. * &dep_var. (0) = &main_spline_var._RCS_lin 
					&main_spline_var._RCS_S1-&main_spline_var._RCS_S%eval(&nb_knots_msv. - 2)
					%do i = 1 %to &nb_osv.;
					%let osv_curr = &&oth_spline_var&i.;
					&osv_curr._RCS_lin &osv_curr._RCS_S1-&osv_curr._RCS_S%eval(&&nb_knots_osv&i. - 2) 
					%end;
					&adjust_var. / covb;

		output out = &outfile_for_analysis. survival = survival;

		%if &nb_knots_msv. = 3 %then %do;
		Overall_association: TEST  &main_spline_var._RCS_lin,&main_spline_var._RCS_S1;
		Non_lin_association: TEST  &main_spline_var._RCS_S1;
		%end;

		%if &nb_knots_msv. = 4 %then %do;
		Overall_association: TEST  &main_spline_var._RCS_lin,&main_spline_var._RCS_S1,&main_spline_var._RCS_S2;
		Non_lin_association: TEST  &main_spline_var._RCS_S1,&main_spline_var._RCS_S2;
		%end;

		%if &nb_knots_msv. = 5 %then %do;
		Overall_association: TEST  &main_spline_var._RCS_lin,&main_spline_var._RCS_S1,&main_spline_var._RCS_S2,&main_spline_var._RCS_S3;
		Non_lin_association: TEST  &main_spline_var._RCS_S1,&main_spline_var._RCS_S2,&main_spline_var._RCS_S3;
		%end;

	&prgm_statements_cox.;

	run;

	quit;
	ods output close;

 	%end;

	/* ----------------- */
	/* Linear regression */
	/* ----------------- */

	%if &typ_reg. = lin %then %do;

	%let text = "Difference in &dep_var. where the reference for &main_spline_var. is &ref_val_FA.";
	%let Yaxis_gplot = "Diff in &dep_var. where the ref value for &main_spline_var. is &ref_val_FA.";

	title1 " ";
	%if %length(%str(&by_factor.)) = 0 %then %do;title2 "   Linear regression model including '&main_spline_var.' using a RCS function with &nb_knots_msv. knots";%end;
	%if %length(%str(&by_factor.)) > 0 %then %do;title2 "   Linear regression model including '&main_spline_var.' using a RCS function with &nb_knots_msv. knots for &by_factor. = &&val_by_factor_&no_cl_byfact.";%end;
	title3 " ";
	title4 "             * - * - * - * - * - * - * - * - * - * - * - * - *";
	title5 " ";

	 	%if %length(%str(&subject_var.)) = 0 %then %do;
		ods output Modelfit = data_fit;
		ods output ParameterEstimates = data_parameters;
		ods output CovB = data_mat_cov;
		%if &print_covar_mat. = 0 %then %do;ods exclude CovB;%end;
		%end;

	 	%if %length(%str(&subject_var.)) > 0 %then %do;
		ods output GEEEmpPEst = data_parameters;
		ods output GEERCov =    data_mat_cov;
		%if &print_covar_mat. = 0 %then %do;ods exclude GEENCov CovB GEERCov;%end;
		%end;


	proc genmod data = &outfile_for_analysis. ;
	%if %length(%str(&subject_var.)) > 0 %then %do;class &subject_var.;%end;
	model &dep_var. = 	&main_spline_var._RCS_lin &main_spline_var._RCS_S1-&main_spline_var._RCS_S%eval(&nb_knots_msv. - 2) 
						%do i = 1 %to &nb_osv.;
						%let osv_curr = &&oth_spline_var&i.;
						&osv_curr._RCS_lin &osv_curr._RCS_S1-&osv_curr._RCS_S%eval(&&nb_knots_osv&i. - 2) 
						%end;

						&adjust_var. / dist = normal link = identity covb;

		output out = &outfile_for_analysis. predicted = Pred_&dep_var. lower = Lo_Pred_&dep_var. upper = Hi_Pred_&dep_var.;

	 	%if %length(%str(&subject_var.)) > 0 %then %do;
		repeated subject = &subject_var. / type = &work_corr_matrix. CovB;
		%end;

		%if &nb_knots_msv. = 3 %then %do;
		contrast "Overall_association" &main_spline_var._RCS_lin 1,&main_spline_var._RCS_S1 1 / wald;
		contrast "Non_lin_association" &main_spline_var._RCS_S1 1 / wald;
		%end;

		%if &nb_knots_msv. = 4 %then %do;
		contrast "Overall_association" &main_spline_var._RCS_lin 1,&main_spline_var._RCS_S1 1,&main_spline_var._RCS_S2 1 / wald;
		contrast "Non_lin_association" &main_spline_var._RCS_S1 1, &main_spline_var._RCS_S2 1 / wald;
		%end;

		%if &nb_knots_msv. = 5 %then %do;
		contrast "Overall_association" &main_spline_var._RCS_lin 1,&main_spline_var._RCS_S1 1,&main_spline_var._RCS_S2 1,&main_spline_var._RCS_S3 1 / wald;
		contrast "Non_lin_association" &main_spline_var._RCS_S1 1,&main_spline_var._RCS_S2 1,&main_spline_var._RCS_S3 1 / wald;
		%end;

	run;

	quit;
	ods output close;

	%end;

	/* ---------------------------------------------------- */
	/* Preparation of data confidence interval calculations */
	/* ---------------------------------------------------- */

	%if &typ_reg. = lin or &typ_reg. = log %then %do;

	data mat_cov (keep = Prm2-Prm&nb_knots_msv.);
	set Data_mat_cov;
	if 1 < _n_ <= &nb_knots_msv.;
	run;

	data _null_;
	set data_parameters;
	if _n_ = 2 then call symput ("P_var_lin",estimate);
	run;

 		%do i = 1 %to %eval(&nb_knots_msv. - 2);
 		data _null_;
 		set data_parameters;
 		if _n_ = (&i. + 2) then call symput ("P_var_S&i.",estimate);
 		run;
 		%end;

	%end;

	%if &typ_reg. = cox %then %do;	

	data mat_cov (keep = Prm2-Prm&nb_knots_msv.);
	set Data_mat_cov;
	if 1 <= _n_ < &nb_knots_msv.;
	rename &main_spline_var._RCS_lin = Prm2;
		%do i = 3 %to &nb_knots_msv.;
		rename &main_spline_var._RCS_S%eval(&i. - 2) = Prm&i.;
		%end;
	run;

	data _null_;
	set data_parameters;
	if _n_ = 1 then call symput ("P_var_lin",estimate);
	run;

 		%do i = 1 %to %eval(&nb_knots_msv. - 2);
 		data _null_;
 		set data_parameters;
 		if _n_ = (&i. + 1) then call symput ("P_var_S&i.",estimate);
 		run;
 		%end;

	%end;


data _null_;
set data_parameters end = final;
retain pb_regression 0;
if &test_name. = . then pb_regression = 1;
if final then do;
if pb_regression = 0 then call symput("Pb_regression",0);
if pb_regression = 1 then call symput("Pb_regression",1);
end;
run;


/* -------------------- */
/* Coordinates of knots */
/* -------------------- */


	%do i = 1 %to &nb_knots_msv.;
	%let V_knot&i. = &&V_knot&i._msv;
	%end;

	%if &pb_regression. = 0 %then %do;	/* Loop with no pb with the regression */

	%let nb_knots_minus1 = %eval(&nb_knots_msv. - 1);

	data export_results;
	do i = 1 to 5;

		%do i_bis = 1 %to &nb_knots_msv.;

	 	if i = &i_bis. then do;
	 	valeur = &&V_knot&i_bis.;

		 S0 = &&V_knot&i_bis. - &ref_val_FA.;

			%do i = 1 %to %eval(&nb_knots_msv. - 2);
			S&i. =	(&&V_knot&i_bis. > &&V_knot&i.)*(&&V_knot&i_bis. - &&V_knot&i.)**3 -
					(&&V_knot&i_bis. > &&V_knot&nb_knots_minus1.) * ((&&V_knot&nb_knots_msv. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&&V_knot&i_bis. - &&V_knot&nb_knots_minus1.)**3 +
					(&&V_knot&i_bis. > &&V_knot&nb_knots_msv.) * ((&&V_knot&nb_knots_minus1. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&&V_knot&i_bis. - &&V_knot&nb_knots_msv.)**3 -
				(
					(&ref_val_FA. > &&V_knot&i.)*(&ref_val_FA. - &&V_knot&i.)**3 -
					(&ref_val_FA. > &&V_knot&nb_knots_minus1.) * ((&&V_knot&nb_knots_msv. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&ref_val_FA. - &&V_knot&nb_knots_minus1.)**3 +
					(&ref_val_FA. > &&V_knot&nb_knots_msv.) * ((&&V_knot&nb_knots_minus1. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&ref_val_FA. - &&V_knot&nb_knots_msv.)**3
				)
					;
			%end;

			%if &exp_beta. = 0 %then %do;
			%if &nb_knots_msv. = 3 %then %do;Y_coordinate = &P_var_lin.*S0 + &P_var_s1.*S1;%end;
			%if &nb_knots_msv. = 4 %then %do;Y_coordinate = &P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2;%end;
			%if &nb_knots_msv. = 5 %then %do;Y_coordinate = &P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2 + &P_var_s3.*S3;%end;
			%end;

			%if &exp_beta. = 1 %then %do;
			%if &nb_knots_msv. = 3 %then %do;Y_coordinate = exp(&P_var_lin.*S0 + &P_var_s1.*S1);%end;
			%if &nb_knots_msv. = 4 %then %do;Y_coordinate = exp(&P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2);%end;
			%if &nb_knots_msv. = 5 %then %do;Y_coordinate = exp(&P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2 + &P_var_s3.*S3);%end;
			%end;

	 	end;

		%end;

	output;
	end;
	run;

	data export_results (drop = i S0 S1);
	set export_results;
		%if &nb_knots_msv. = 3 %then %do;
			if _n_ in (4,5) then do;
			valeur = .;
			Y_coordinate = .;
			end;
		%end;
		%if &nb_knots_msv. = 4 %then %do;
		drop S2;
			if _n_ = 5 then do;
			valeur = .;
			Y_coordinate = .;
			end;
		%end;
	%if &nb_knots_msv. = 5 %then %do;drop S2 S3;%end;
	run;

	data _null_;
	inc = (&max_Xaxis_FA. - &min_Xaxis_FA.)/200;
	call symput ("inc",inc);
	run;


/* ---------------------------------------- */
/* Dose-response curve calculation + 95% CI */
/* ---------------------------------------- */

	data temp_estimations1;
	do varX = &min_Xaxis_FA. to &max_Xaxis_FA. by &inc.;
	S0 = varX - &ref_val_FA.;

		%do i = 1 %to %eval(&nb_knots_msv. - 2);
		S&i. =		(varX > &&V_knot&i.)*(varX - &&V_knot&i.)**3 -
					(varX > &&V_knot&nb_knots_minus1.) * ((&&V_knot&nb_knots_msv. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (varX - &&V_knot&nb_knots_minus1.)**3 +
					(varX > &&V_knot&nb_knots_msv.) * ((&&V_knot&nb_knots_minus1. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (varX - &&V_knot&nb_knots_msv.)**3 -
				(
					(&ref_val_FA. > &&V_knot&i.)*(&ref_val_FA. - &&V_knot&i.)**3 -
					(&ref_val_FA. > &&V_knot&nb_knots_minus1.) * ((&&V_knot&nb_knots_msv. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&ref_val_FA. - &&V_knot&nb_knots_minus1.)**3 +
					(&ref_val_FA. > &&V_knot&nb_knots_msv.) * ((&&V_knot&nb_knots_minus1. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&ref_val_FA. - &&V_knot&nb_knots_msv.)**3
				)
					;
		%end;

	%if &nb_knots_msv. = 3 %then %do;varY = &P_var_lin.*S0 + &P_var_s1.*S1;output;%end;
	%if &nb_knots_msv. = 4 %then %do;varY = &P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2;output;%end;
	%if &nb_knots_msv. = 5 %then %do;varY = &P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2 + &P_var_s3.*S3;output;%end;

	end;
	run;

	proc iml;
		%do i = 1 %to 200;

			if &nb_knots_msv. = 3 then do;
			use temp_estimations1;read point &i. var {S0 S1} into mat_x;close temp_estimations1;
			use mat_cov;read all var {Prm2 Prm3} INTO matcov;close mat_cov;
			end;
			if &nb_knots_msv. = 4 then do;
			use temp_estimations1;read point &i. var {S0 S1 S2} into mat_x;close temp_estimations1;
			use mat_cov;read all var {Prm2 Prm3 Prm4} INTO matcov;close mat_cov;
			end;
			if &nb_knots_msv. = 5 then do;
			use temp_estimations1;read point &i. var {S0 S1 S2 S3} into mat_x;close temp_estimations1;
			use mat_cov;read all var {Prm2 Prm3 Prm4 Prm5} INTO matcov;close mat_cov;
			end;

		t_mat_x = T(mat_x);
		produit = mat_x*matcov*t_mat_x;

		%if &i. = 1 %then %do;mat_final=produit;%end;
		%if &i. > 1 %then %do;mat_final=mat_final//produit;%end;

		create Data_product_mat from mat_final;append from mat_final;close Data_product_mat;

		%end;

	quit;

	data estimations (drop = S0-S%eval(&nb_knots_msv. - 2));
	merge Data_product_mat temp_estimations1;
	run;

	data estimations (drop = produit col1 typ_reg);
	set estimations;
	produit = col1;
	typ_reg = "%sysfunc(lowcase(%str(&typ_reg.)))";
		if typ_reg = "lin" then do;
		Lower_CL = varY - 1.96*sqrt(produit);
		Upper_CL = varY + 1.96*sqrt(produit);
		end;
		if typ_reg in ("log","cox") then do;
			if &exp_beta. = 0 then do;
			Lower_CL = varY - 1.96*sqrt(produit);
			Upper_CL = varY + 1.96*sqrt(produit);
			end;
			if &exp_beta. = 1 then do;
			varY = exp(varY);
			Lower_CL = exp(log(varY) - 1.96*sqrt(produit));
			Upper_CL = exp(log(varY) + 1.96*sqrt(produit));
			end;
		end;
	rename varY = Estimation;
	run;


	/* ------------------------------------------------------------- */
	/* Estimation of the association for a specific value (+ 95% CI) */
	/* ------------------------------------------------------------- */

		%if %length(%str(&specif_val.)) > 0 %then %do;

			%do sp_values_index = 1 %to &nb_sp_values.;
			%let sp_values&sp_values_index. = %scan(&specif_val.,&sp_values_index.,%STR( ));


			data temp_est_valpart;
			varX = &&sp_values&sp_values_index.;
			S0 = varX - &ref_val_FA.;

				%do i = 1 %to %eval(&nb_knots_msv. - 2);
				S&i. =	(varX > &&V_knot&i.)*(varX - &&V_knot&i.)**3 -
						(varX > &&V_knot&nb_knots_minus1.) * ((&&V_knot&nb_knots_msv. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (varX - &&V_knot&nb_knots_minus1.)**3 +
						(varX > &&V_knot&nb_knots_msv.) * ((&&V_knot&nb_knots_minus1. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (varX - &&V_knot&nb_knots_msv.)**3 -
					(
						(&ref_val_FA. > &&V_knot&i.)*(&ref_val_FA. - &&V_knot&i.)**3 -
						(&ref_val_FA. > &&V_knot&nb_knots_minus1.) * ((&&V_knot&nb_knots_msv. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&ref_val_FA. - &&V_knot&nb_knots_minus1.)**3 +
						(&ref_val_FA. > &&V_knot&nb_knots_msv.) * ((&&V_knot&nb_knots_minus1. - &&V_knot&i.)/(&&V_knot&nb_knots_msv. - &&V_knot&nb_knots_minus1.)) * (&ref_val_FA. - &&V_knot&nb_knots_msv.)**3
					)
						;
				%end;

			%if &nb_knots_msv. = 3 %then %do;varY = &P_var_lin.*S0 + &P_var_s1.*S1;output;%end;
	 		%if &nb_knots_msv. = 4 %then %do;varY = &P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2;output;%end;
		 	%if &nb_knots_msv. = 5 %then %do;varY = &P_var_lin.*S0 + &P_var_s1.*S1 + &P_var_s2.*S2 + &P_var_s3.*S3;output;%end;

			run;

			proc iml;

				if &nb_knots_msv. = 3 then do;
				use temp_est_valpart;read point 1 var {S0 S1} into mat_x;close temp_est_valpart;
				use mat_cov;read all var {Prm2 Prm3} INTO matcov;close mat_cov;
				end;
				if &nb_knots_msv. = 4 then do;
				use temp_est_valpart;read point 1 var {S0 S1 S2} into mat_x;close temp_est_valpart;
				use mat_cov;read all var {Prm2 Prm3 Prm4} INTO matcov;close mat_cov;
				end;
				if &nb_knots_msv. = 5 then do;
				use temp_est_valpart;read point 1 var {S0 S1 S2 S3} into mat_x;close temp_est_valpart;
				use mat_cov;read all var {Prm2 Prm3 Prm4 Prm5} INTO matcov;close mat_cov;
				end;

			t_mat_x = T(mat_x);
			produit = mat_x*matcov*t_mat_x;
			create var_est_valpart var {produit};append;close var_est_valpart;
			run;
			quit;

			data est_valpart (drop = S0-S%eval(&nb_knots_msv. - 2));
			merge var_est_valpart temp_est_valpart;
			run;

			data data_sp_values (drop = produit typ_reg);
			set est_valpart;
			typ_reg = "%sysfunc(lowcase(%str(&typ_reg.)))";
	 			if typ_reg = "lin" then do;
				Lower_CL = varY - 1.96*sqrt(produit);
				Upper_CL = varY + 1.96*sqrt(produit);
				end;
				if typ_reg in ("log","cox") then do;
					if &exp_beta. = 0 then do;
					Lower_CL = varY - 1.96*sqrt(produit);
					Upper_CL = varY + 1.96*sqrt(produit);
					end;
					if &exp_beta. = 1 then do;
					varY = exp(varY);
					Lower_CL = exp(log(varY) - 1.96*sqrt(produit));
					Upper_CL = exp(log(varY) + 1.96*sqrt(produit));
					end;
				end;
				%if %length(%str(&round.)) > 0 %then %do;
				varY = round(varY,&round.);
				Lower_CL = round(Lower_CL,&round.);
				Upper_CL = round(Upper_CL,&round.);
				%end;
			run;

			data data_sp_values;
			set data_sp_values;
			rename varX = &main_spline_var.;
			rename varY = &for_rename.;
			run;
			
				%if %length(%str(&by_factor.)) = 0 %then %do;

					%if &sp_values_index. = 1 %then %do;
					data list_specif_values;
					set data_sp_values;
					run;
					%end; 

					%if &sp_values_index. > 1 %then %do;
					data list_specif_values;
					set list_specif_values data_sp_values;
					run;
					%end; 

				%end;

				%if %length(%str(&by_factor.)) > 0 %then %do;

					%if &sp_values_index. = 1 %then %do;
					data list_specif_values_&no_cl_byfact.;
					set data_sp_values;
					run;
					%end; 

					%if &sp_values_index. > 1 %then %do;
					data list_specif_values_&no_cl_byfact.;
					set list_specif_values_&no_cl_byfact. data_sp_values;
					run;
					%end; 

				%end;

			%end;

		%end;


	/* ------------------ */
	/* Plot of the curves */
	/* ------------------ */

		%if &no_graph. = 0 %then %do;

		data for_gplot1 (keep = valeur knots);
		set export_results;
		Knots = Y_coordinate;
		run;

		data for_gplot2;
		set estimations for_gplot1;
		run;

		data for_gplot3;
		set for_gplot2;
		if (valeur = . or (&min_Xaxis_FA. <= valeur <= &max_Xaxis_FA. ));
		run;

		proc means data = for_gplot3 noprint min max;
		output out = for_plot_obs	min (lower_cl) = mini_low_cl max (upper_cl) = maxi_upp_cl;
		run;

		data _null_;
		set for_plot_obs;
		call symput("min_Y",mini_low_cl);
		call symput("max_Y",maxi_upp_cl);
		run;

			%if (&typ_reg. = lin or &typ_reg. = log) %then %do;
			data for_gplot4 (drop = Pred_&dep_var.);
			set for_gplot3 &outfile_for_analysis. 
				(	keep = &main_spline_var._RCS_lin Pred_&dep_var.
					where = ((&min_Xaxis_FA. <= &main_spline_var._RCS_lin <= &max_Xaxis_FA.)
								and (Pred_&dep_var. ne .)));
			Obs = &min_Y. - (&max_Y. - &min_Y.)/50 - (&width_band_obs. * (rand("uniform") * (&max_Y. - &min_Y.)/100));
			run;
			%end;

			%if &typ_reg. = cox %then %do;
			data for_gplot4 (drop = survival);
			set for_gplot3 &outfile_for_analysis. 
				(	keep = &main_spline_var._RCS_lin survival
					where = ((&min_Xaxis_FA. <= &main_spline_var._RCS_lin <= &max_Xaxis_FA.)
								and (survival ne .)));
			Obs = &min_Y. - (&max_Y. - &min_Y.)/50 - (&width_band_obs. * (rand("uniform") * (&max_Y. - &min_Y.)/100));
			run;
			%end;

		%let legend = ;
		%if &no_legend. = 1 %then %let legend = noautolegend;

		proc sgplot data = for_gplot4 &legend.;
		%if &no_title. = 1 %then %do;title1 " ";%end;
			%if &no_title. = 0 %then %do;
				%if %length(%str(&by_factor.)) = 0 %then %do;
				title1 &titre_gplot.;
				%end;
				%if %length(%str(&by_factor.)) > 0 %then %do;
				title1 &titre_gplot.;
				title2 &titre_gplot_BF.;
				%end;
			%end;
		series X=varX Y=Estimation / LINEATTRS=(color = red THICKNESS=2);
		series X=varX Y=lower_CL / LINEATTRS=(color = black PATTERN=3);
		series X=varX Y=upper_CL / LINEATTRS=(color = black PATTERN=3);
		%if &display_knots. = 1 %then %do;scatter X=valeur  Y=knots / MARKERATTRS=(color=red size=8 SYMBOL=circlefilled);%end;
		%if &print_obs. = 1 %then %do;scatter X=&main_spline_var._RCS_lin Y=obs / MARKERATTRS=(color=black size=4 SYMBOL=plus);%end;
		%if &X_ref_line. = 1 %then %do;refline &ref_val_FA. / axis = X LINEATTRS = (color = green PATTERN=3);%end;
			%if &Y_ref_line. = 1 %then %do;
				%if &exp_beta. = 0 %then %do;refline 0 / axis = Y LINEATTRS = (color = green PATTERN=3);%end;
				%if &exp_beta. = 1 %then %do;refline 1 / axis = Y LINEATTRS = (color = green PATTERN=3);%end;
			%end;
		%if &no_label_X. = 0 %then %do;xaxis label = "&main_spline_var.";%end;
		%if &no_label_X. = 1 %then %do;xaxis display = (nolabel);%end;
		%if &no_label_Y. = 0 %then %do;yaxis label = &Yaxis_gplot.;%end;
		%if &no_label_Y. = 1 %then %do;yaxis display = (nolabel);%end;
		run;

		%end;

	/* ------------------------------- */
	/* Export to txt format for curves */
	/* ------------------------------- */

		%if &export_curves. = 1 %then %do;

			%if %length(%str(&by_factor.)) = 0 %then %do;

				%if &print_obs. = 0 %then %do;
				data export_curves (keep = &main_spline_var. estimation Lower_CL Upper_CL 
								where = (&main_spline_var. ne .));
				set for_gplot4;
				rename varX = &main_spline_var.;
				run;
				%end;

				%if &print_obs. = 1 %then %do;
				data export_curves (keep = &main_spline_var. estimation Lower_CL Upper_CL Obs_X Obs_Y
								where = (&main_spline_var. ne . or Obs_X ne .));
				set for_gplot4;
				if &main_spline_var._RCS_lin = . then Obs = .;
				rename varX = &main_spline_var.;
				rename &main_spline_var._RCS_lin = Obs_X;
				rename Obs = Obs_Y;
				run;
				%end;

			PROC EXPORT DATA= export_curves
	        OUTFILE= "&dir_export_curves.\export_curves_rcs.txt" 
	        DBMS=TAB REPLACE;
	     	PUTNAMES=YES;
			RUN;

			%end;

			%if %length(%str(&by_factor.)) > 0 %then %do;

				%if &print_obs. = 0 %then %do;
				data export_curves_&no_cl_byfact. (	keep = &main_spline_var. estimation Lower_CL Upper_CL 
												where = (&main_spline_var. ne .));
				set for_gplot4;
				rename varX = &main_spline_var.;
				run;
				%end;

				%if &print_obs. = 1 %then %do;
				data export_curves_&no_cl_byfact. (	keep = &main_spline_var. estimation Lower_CL Upper_CL Obs_X Obs_Y
												where = (&main_spline_var. ne . or Obs_X ne .));
				set for_gplot4;
				if &main_spline_var._RCS_lin = . then Obs = .;
				rename varX = &main_spline_var.;
				rename &main_spline_var._RCS_lin = Obs_X;
				rename Obs = Obs_Y;
				run;
				%end;

			PROC EXPORT DATA= export_curves_&no_cl_byfact.
	        OUTFILE= "&dir_export_curves.\export_curves_rcs_&no_cl_byfact..txt" 
	        DBMS=TAB REPLACE;
	     	PUTNAMES=YES;
			RUN;

			%end;

		%end;


	/* --------- */
	/* AIC value */
	/* --------- */

	%if (&typ_reg. = lin or &typ_reg. = log) and %length(%str(&subject_var.)) = 0 %then %do;

	data _null_;
	set data_fit;
	if criterion = "AIC (smaller is better)" then call symput("AIC",value);
	run;

	%end;

	%if &typ_reg. = cox %then %do;

	data _null_;
	set data_fit;
	if criterion = "AIC" then call symput("AIC",WithCovariates);
	run;

	%end;

	/* ----------------------------------- */
	/* ORs or HRs in the SAS OUTPUT window */
	/* ----------------------------------- */

		%if &exp_beta. = 1 %then %do;

			%if &print_OR_HR. = 1 %then %do;

				%if %length(%str(&by_factor.)) = 0 %then %do;

				data List_OR_HR (drop = valeur knots);
				set for_gplot2;
				if varX ne .;
				rename varX = &main_spline_var.;
				rename estimation = &for_rename.;
				run;

				Proc print data = List_OR_HR noobs;
				title1 " ";
				title2 "------------------------------------------------------------------------------------------------------------- ";
				title3 " ";
				title4 "List of the &for_rename [95% CI] for values of &main_spline_var. ranging from &min_Xaxis_FA. to &max_Xaxis_FA. (ref value for &main_spline_var.: &ref_val_FA.)";
				title5 " ";
				run; 

				%end;	

				%if %length(%str(&by_factor.)) > 0 %then %do;

				data List_OR_HR_&no_cl_byfact. (drop = valeur knots);
				set for_gplot2;
				if varX ne .;
				rename varX = &main_spline_var.;
				rename estimation = &for_rename.;
				run;

				Proc print data = List_OR_HR_&no_cl_byfact. noobs;
				title1 " ";
				title2 "------------------------------------------------------------------------------------------------------------- ";
				title3 " ";
				title4 "List of the &for_rename [95% CI] for values of &main_spline_var. ranging from &min_Xaxis_FA. to &max_Xaxis_FA. (ref value for &main_spline_var.: &ref_val_FA.)";
				title5 "For &by_factor. = &&val_by_factor_&no_cl_byfact.";
				title6 " ";
				run; 

				%end;	

			%end;

		%end;

	/* ------------------------------------------------ */
	/* List of the specific values in the OUTPUT window */
	/* ------------------------------------------------ */

		%if %length(%str(&by_factor.)) = 0 %then %do;

			%if %length(%str(&specif_val.)) > 0 %then %do;
			proc print data = List_specif_values noobs;
			title1 " ";
			title2 "---------------------------------------------------------------------------------------------------- ";
			title3 " ";
			title4 "List of the &for_rename [95% CI] for specific values of &main_spline_var. (ref value for &main_spline_var.: &ref_val_FA.)";
			title5 " ";
			run;
			%end;

		%end;

		%if %length(%str(&by_factor.)) > 0 %then %do;

			%if %length(%str(&specif_val.)) > 0 %then %do;
			proc print data = list_specif_values_&no_cl_byfact. noobs;
			title1 " ";
			title2 "---------------------------------------------------------------------------------------------------- ";
			title3 " ";
			title4 "List of the &for_rename [95% CI] for specific values of &main_spline_var. (ref value for &main_spline_var.: &ref_val_FA.)";
			title5 "For &by_factor. = &&val_by_factor_&no_cl_byfact.";
			title6 " ";
			run;
			%end;

		%end;

	%end;	/* End of loop with no pb with the regression */

%end;

/* ---------------------- */
/* Delete temporary files */
/* ---------------------- */

	%if &no_delete_files. = 0 %then %do;
	proc datasets library = work nolist;
	delete 	export_results temp mat_cov estimations temp_estimations1 for_gplot1 for_gplot2 for_gplot3 for_gplot4
				Data_sp_values
				temp_est_valpart Est_valpart Var_est_valpart Data_mat_cov data_parameters Data_product_mat
				&infile._for_analysis data_fit for_plot_obs;
	run;
	%end;

	%if %length(%str(&by_factor.)) > 0 %then %do;
	proc datasets library = work nolist;
	delete Data_by_factor;
	run;
	%end;

proc datasets library = work nolist;	 
delete  out_univariate out_distrib Defknots_from_univ ; 
run;


/* ---------------------------------- */
/* Misc information in the LOG Window */
/* ---------------------------------- */

%put ;
%put &message_add1. ;
%put &message_add2. ;
%put ;
%put Spline variables have been created in the "&outfile_for_analysis." SAS datafile, located in the SAS Working library;
%if &no_graph. = 1 %then %put No graph has been created;
%if %length(%str(&message_ref_val.)) > 0 %then %put &message_ref_val.;
%put ;
%put - - - - - - -;
%put ;
%put Name of the %eval(&nb_knots_msv. - 1) spline variables of &main_spline_var. that have just been created:;
%put ;
%put Spline #1 = &&main_spline_var._RCS_lin;
%do i = 1 %to %eval(&nb_knots_msv. - 2);
%put Spline #%eval(&i. + 1) = &&main_spline_var._RCS_S&i.;
%end;
%put ;
%put - - - - - - -;
%do i = 1 %to &nb_osv.;
%let name_osv_curr = &&oth_spline_var&i.;
%put ;
%put Name of the %eval(&&nb_knots_osv&i. - 1) spline variables of &name_osv_curr. that have just been created:;
%put ;
%put Spline #1 = &name_osv_curr._RCS_lin;
	%do j = 1 %to %eval(&&nb_knots_osv&i. - 2);
	%put Spline #%eval(&j. + 1) = &&name_osv_curr._RCS_S&j.;
	%end;
%put ;
%put - - - - - - -;
%end;
%put ;

%if %length(%str(&typ_reg.)) > 0 %then %do;
	%if %length(%str(&subject_var.)) = 0 %then %do;
	%put The AIC value for the model is: %trim(%left(&AIC.));
	%end;
%end;


%put ;
%put ;
%put ------------------------- End of miscellaneous information ------------------------;
%put ;
%put ;
%put ;
%put ;

%if &pb_regression. = 1 %then %do;

%put ;
%put ;
%put *		+------------------+;
%put *		| **** ERROR! **** |;
%put *		+------------------+;
%put ;
%put The regression has encountered some problems. Please check (for instance, the range of &main_spline_var. may be too wide compared to the 5th or the 95th percentile).;
%put ;
%put ;

%end;


/* End of PUT instructions in the LOG Window */


%end;	/* End of mistake-free loop #2 */

%end;	/* End of BY_FACTOR loop */

proc datasets library = work nolist;
delete 	&infile._II;
run;

option notes;
title ;
title1 ;
title2 ;
title3 ;
title4 ;
title5 ;
title6 ;

quit;

%mend RCS_Reg;

/* Counting the number of words in a string */
%MACRO WORDS(STRING);
%LOCAL COUNT WORD;
%LET COUNT=1;
	%LET WORD=%SCAN(&STRING.,&COUNT.,%STR( ));
	%DO %WHILE(&WORD. NE );
	%LET COUNT=%EVAL(&COUNT. + 1);
	%LET WORD=%SCAN(&STRING.,&COUNT.,%STR( ));
	%END;
%EVAL(&COUNT. - 1)
%MEND WORDS;

