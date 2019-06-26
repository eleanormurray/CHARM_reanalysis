/*****************************************************************************************
* ADHERENCE ADJUSTMENT IN PLACEBO-CONTROLLED RANODMIZED TRIALS:							 *
* AN APPLICATION TO THE CANDESARTAN IN HEART FAILURE RANDOMIZED TRIAL 					 *
* Authors: Eleanor Murray (ejmurray@bu.edu), Brian Claggett, Bradi Granger,				 *
*          Scott Solomon, Miguel Hernan                                                  *
*                                                                                        *
*                                                                                        *
* Version May 2019	.						                                             *    
*                                                                                        *
******************************************************************************************/


/*
Copyright (c) 2014, 2018, The President and Fellows of Harvard College
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.
This software is provided under the standard MIT License:
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
*/

libname charm "<long-format data>";
libname master "<wide-format data>";

/*******************************************************/
/************Set up common macro variables**************/
/*******************************************************/

%let eof = 182; /*3.5 years*/

%let covs1 = age0 ef0 comorbs_count0 nyha_b hr0 sex meds_count0 dbp0 sbp0 bminum0 curr_smoke0 hxpaim0; /*Original Granger et al list*/

%let covs2 = age0 ef0  nyha_b3 nyha_b4 hr0 sex study_ dbp0 sbp0 bminum0 curr_smoke0 	
		dgdiu_0 dgbb_0 dglild_0 dgacea_0 dgacei_0 dgdig_0 dgvas_0 dgantiarr_0 dgoral_0 dgantip_0  dgothercv_0
		hxaf_0 hxap_0 hxcabg_0 hxcan_0 hxdiam_0 hxhyp_0 hxicd_0 hxpcor_0 hxstr_0 hxpmi_0 hxpaim_0; /*Expanded list with dummy variables for categorical covariates*/

%let covs_t = dgdig dgdiu dgbb dgccb dgvas dgantiarr dglild dgoral dgacea dgantip dgacei dgothercv
				cabgh diab icdh pcorh strh paimh curr_smoke dbp_t sbp_t hrate_t nyha2 nyha3 nyha4; /*Post-randomization covariates*/
				
proc format;
value adher 
	1 = "A: >80% adherent"
	0 = "B: =<80% adherent"
	;
value death 
	1 = "A: died"
	0 = "B: alive"
;
value rand
	1 = "Candesartan"
	0 = "placebo"
;
run;

/*******************************************************/
/*Macro cumavg: Create adherence & censoring covariates*/
/*******************************************************/
%macro cumavg(inset = , outset = ,  varname_out = , varout_baseline = , id = , time = ,	 completecase = 0, compyes = 0, 
		n_cens_cond = , censvarname = , censif_var1 =,  censif_cond1 = ,  censif_var2=, censif_cond2 = , nocensbefore = ,
		recode = , recode_var1 = , recodeifeq_cond1 = , recodeto_cond1 = ,  recodeifeq_cond2 = , recodeto_cond2 = ,  recodeifeq_cond3 = , recodeto_cond3 = ,  norecodebefore = ,
	 	LTFUvarname = , LTFUaftervisit = , LTFU = 0, LTFUcond_var1=  , dichot = 0.8, outdatatype = 0, restrict_arms = 0);

proc sort data = &inset;
by &id &time;

%if &completecase = 1 %then %do;

	%if &compyes = 1 %then %do;
		%let varname = compyes;
	%end;
	%else %if &compyes = 0 %then %do;
		%let varname = compl;
	%end;

	data events;
	set &inset;
		if visit_t = 1 then output;
	run;

	proc means data = events mean noprint;
	by &id;
		var &varname;
		output out = temp (keep = &id cumavg_eof_cc) mean = cumavg_eof_cc;
	run;

	%if &dichot ne . %then %do;	
		%let binavg = &varname_out._bin;
	%end;

	data &outset;
	merge &inset temp;
	by &id;

		&varname_out = cumavg_eof_cc;
	
		if &dichot ne . then do;	
		 	if &varname_out > &dichot then &binavg = 1; else if . < &varname_out le &dichot then &binavg = 0; 
		end;
	run;

	
%end;
%else %do;
	%if &compyes = 0 %then %do;
		%let var_t = compl_t;
		%let var_t0 = compl_t0;
	%end;
	%else %if &compyes = 1 %then %do;
		%let var_t = compyes_t;
		%let var_t0 = compyes_t0;
	%end;

	data &outset;
	set &inset;
	by &id &time;

		retain &varout_baseline;
		%if &n_cens_cond > 0 %then %do;
			retain &censvarname;
		%end;	
		%if &LTFU =1 %then %do;
			retain &LTFUvarname;
		%end;

		if first.b_patid then do;

			&varout_baseline = &&var_t0;
			&varname_out 	 = &&var_t;

			/*create censoring not based on loss to follow-up*/
			%if &n_cens_cond = 1  %then %do; 
				if &censif_var1 = &censif_cond1 	then &censvarname = 1;
				else &censvarname = 0;
			%end;
			%else %if &n_cens_cond =2 %then %do; 
				if &censif_var1 = &censif_cond1 	then &censvarname = 1;
				else if &censif_var2 =  &censif_cond2 	then &censvarname = 1;
				else &censvarname = 0;	
			%end;
			/*code as adherent if discontinued due to allowed reasons*/
			
			%if &recode = 1 %then %do;
				if &recode_var1 = &&recodeifeq_cond1 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond1;
				end;
				else &varname_out = &&var_t;
			%end;
			%if &recode = 2 %then %do;
				if &recode_var1 = &&recodeifeq_cond1 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond1;
					recode_flag_t = 1;
				end;
				else if &recode_var1 = &&recodeifeq_cond2 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond2;
					recode_flag_t = 1;
				end;
				else if &recode_var1 = &&recodeifeq_cond3 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond3;
					recode_flag_t = 1;
				end;
				else &varname_out = &&var_t;
			%end;
		
			/*create censoring based on loss to follow-up*/
			%if &LTFU = 1 %then %do;
				&LTFUvarname = 0;
			%end;
		end;
		else do;

			&varname_out = &&var_t;
			%if &LTFU = 1 %then %do;
				if comp_m = 1 then do;
					&varname_out = &&var_t;
					if &LTFUvarname ne . then &LTFUvarname = 0;
				end;
				else if comp_m = 0 then do;
				 	if . < &LTFUcond_var1 < &&LTFUaftervisit then do; 
						if &LTFUvarname ne . then &LTFUvarname = 0;
						&varname_out = &&var_t;
					end;
					else if &LTFUcond_var1 = &&LTFUaftervisit   then do;
						&varname_out = &&var_t;
						if &LTFUvarname ne . then &LTFUvarname = 1;
					end;
					else if &LTFUcond_var1 > &&LTFUaftervisit then do;
						&varname_out = .;
						&LTFUvarname = .;
					end;						
				end;
			%end;

			%if &recode = 1 %then %do;
				if &recode_var1 = &&recodeifeq_cond1 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond1;
				end;
				else &varname_out = &&var_t;
			%end;		
			%if &recode = 2 %then %do;
				if &recode_var1 = &&recodeifeq_cond1 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond1;
					recode_flag_t = 1;
				end;
				else if &recode_var1 = &&recodeifeq_cond2 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond2;
					recode_flag_t = 1;
				end;
				else if &recode_var1 = &&recodeifeq_cond3 and weektime > &&norecodebefore then do;
					&varname_out = &recodeto_cond3;
					recode_flag_t = 1;
				end;
				else &varname_out = &&var_t;
			%end;	

			%if &n_cens_cond = 1 %then %do; /*create censoring not based on loss to follow-up*/
				if &censvarname = 1 then do;
					&censvarname = .;
					&varname_out = .;
				end;
				else if &censvarname = . 		then &varname_out = .;
				else if &censif_var1 = &censif_cond1 	then &censvarname = 1;
				else &censvarname = 0;
			%end;
			%else %if &n_cens_cond = 2 %then %do; /*create censoring not based on loss to follow-up*/
				if &censvarname = 1 then do;
					&censvarname = .;
					&varname_out = .;
				end;
				else if &censvarname = . 		then &varname_out = .;
				else if &censif_var1 = &censif_cond1 	then &censvarname = 1;
				else if &censif_var2 = &censif_cond2 	then &censvarname = 1;
				else &censvarname = 0;
			%end;
		end;
	run;

	%let avgname = &varname_out._avg;
	%let avgeof = &avgname._eof;

	data &outset;
	set &outset;
	by &id &time;

		retain &avgname;

		if first.b_patid then do;
			&avgname = &varname_out;
		end;

		else do;
			if visit_t = 0 then do;
				&avgname = &avgname;
			end;
				
			else if visit_t = 1 then do;
				&avgname = (&avgname*newvisit + &varname_out)/(newvisit + 1);
			end;
		end;
	run;

	data last;
	%if &LTFU = 1 %then %do;
		set &outset (where = (&LTFUvarname ne .));
	%end;
	%else %if &LTFU = 0 %then %do;
		set &outset (where = (death ne .));
	%end;
	by &id &time;

		if last.b_patid then do;
			&avgeof = &avgname;
			output;
		end;

		keep b_patid &avgeof;
	run;

	data &outset;
	merge  &outset last;
		by &id;
	run;
%end;

/*output dataset for cumulative incidence analysis*/
%if &outdatatype = 1 %then %do;
	data &outset;

	set &outset (where = (death ne .));
	by &id &time;	
		if last.b_patid then output;
	run;

%end;

/*output dataset with visits only for best replication*/
%else %if &outdatatype = 2 %then %do;
	data &outset;
	set &outset (where = (visit_t = 1));
		by &id &time;
		if last.b_patid then do;
			if death = 1 then event_new = 1;
		end;
	run;
%end;


/*output all data while alive*/
%else %if &outdatatype = 3 %then %do;
	data &outset (where = (event_new ne .));
		set &outset ;
	run;
%end;


/*output all data while alive*/
%else %if &outdatatype = 4 %then %do;
	data &outset;set &outset;
		by &id &time;
		if first.b_patid then do;
			lagevent = 0;
		end;
		else do;
			lagevent = lag(event_new);
		end;
	run;	
	data &outset (where = (event_new ne . & lagevent = 0));
		set &outset ;
	run;
%end;
%mend cumavg;


/**********************************************************************/
/*Macro nostdz: Unadjusted & adjusted analyses without standardization*/
/**********************************************************************/

%macro nostdz(inset = , outdest =, title_main = ,  expvar = , expvar0_bin =, covs= ,  class=, id = b_patid, crude = , recode = 0,
			logit = , start = , stop = 0, event = , itt = , ppe = , ppe_type = , subset = rand , restrict_subset =  ,  censLTFU = , 
			ipcw = , iptw = 0, weights = stabw1, censA = 0, pooled_l = 0, timevar = , ag = 0, boots = 0, nboot = );

title &title_main;

proc sort data=&inset out = replic;
by &subset;
run;

data replic;
set replic;
if &restrict_subset ne . then do;
	if &subset = &restrict_subset then output;
end;

%if &boots = 1 %then %do;
proc sort data = replic;
by b_patid; 
run;
data replic ;
  set replic end = _end_ ;
  by b_patid;
 retain _id ;
  if _n_ = 1 then _id = 0;
  if first.b_patid then do;
  	_id = _id + 1 ;	
  end;
  if _end_ then do ;
     call symput("nids",trim(left(_id)));
  end;
run;
data ids ;
   do bsample = 1 to &nboot;
       do _id = 1 to &nids ;
           output ;
       end;
   end;
run;
proc surveyselect data= ids 
         method = urs
         n= &nids
         seed = 1232  
         out = _idsamples (keep = bsample _id  numberhits  ) 
         outall  noprint  ;       
      strata bsample ;
run;



%end;


data replic;
set replic;

	%if &recode = 1 %then %do;
		if . < &expvar < 2 then adher = 1;
		else if  &expvar ge 2 then adher = 0;
		else if &expvar= . then adher= 0;
	%end;
	%else %if &recode = 2 %then %do;
		adher = 3-&expvar;
	%end;

	%if &recode = 3 %then %do;
		if &expvar > 0.8 then adher = 1;
		else if . <  &expvar =< 0.8 then adher = 0;
		else if  &expvar= . then adher= 0;

	%end;

	%else %if &recode = 4 %then %do;
		if &expvar > 0.8 then adher = 1;
		else if . <  &expvar =< 0.8 then adher = 0;
		if &censLTFU = 1 then adher = .;
	%end;	

run;

%if &censA = 1 %then %do;
	proc sort data = replic; 
		by b_patid &timevar;
	run;
	data replic;
	set replic;
	by b_patid &timevar;
	
	retain censA adher0;
	if first.b_patid then do; 
		censA = 0;
		adher0 = adher;
	end;
	else do;
		if censA = 1 or censA = . then censA = .;
		else if adher0 = adher and censA ne 1 then censA = 0;
		else if adher0 ne adher or censA = 1 then censA = 1;
	end;
	if censA = . then delete;
	run;
	
	proc freq data = replic;
	tables censA / missing list;
	title 'check censA';
	run;
%end;
proc freq data = &inset;
tables &expvar / missing list;
run;
proc freq data = replic;
tables adher &event / missing list;
run;

%if &logit = 1 %then %do;

	proc logistic data=replic descending; 
	%if &subset ne . %then %do; by &subset;%end;
		class sex (ref="1") /param = ref;
		%if &class = 1 %then %do;
			class nyha_b / param = ref;
		%end;
		%if &censA = 1 %then %do; 
			where censA = 0;
		%end;

		%if &crude = 1 %then %do;
			model &event = adher;
		%end;
		%else %if &crude = 0 %then %do;
	
			model &event = adher &covs ;
				units age0 = 10 hr0 = 10 dbp0 = 10 sbp0 = 10;
		%end;
		%if &ipcw = 1 or &iptw = 1 %then %do;	
			weight &weights;
		%end;
			ods output parameterestimates =  param;
			ods output oddsratios = ors;
		run;

	data pvalues;
	set param;
		%if &subset ne . %then %do; 
			keep &subset; 
		%end;
		keep Variable probchisq;

	data ors;
	set ors;
		variable = effect;
		%if &subset ne . %then %do; 
			keep &subset; 
		%end;
		keep variable OddsRatioEst lowercl uppercl;
	%if &subset ne . %then %do;
		proc sort data = ors; by &subset variable;
		proc sort data = pvalues; by &subset variable;
		data all;
		merge ors pvalues;
			by &subset variable;
	%end;
	%else %do;
		proc sort data = ors; by variable;
		proc sort data = pvalues; by variable;
		data all;
		merge ors pvalues;
			by variable;	
	%end;
	proc printto print = &outdest;
	run;

	proc print data = all;
		%if &subset ne . %then %do; 
			by &subset; 
		%end;
		title &title_main;
		%if &subset ne . %then %do; 
			format &subset rand.; 
		%end;		
	run;

	proc printto;
	run;
%end;
%else %if &logit = 0 %then %do;

	proc sort data = replic;
	by b_patid weektime;
	run;

	data replic_surv;
		set replic;
		by b_patid weektime;


		%if &start = 0 %then %do;
			lag_time = lag(weektime);
			if first.b_patid then startweek = .;
			else startweek = lag_time;

			%if &censA = 1 %then %do; 
				where censA = 0;
			%end;

			%if &stop = 0 %then %do;
				stopweek = weektime;
				if &event = 1 then stopweek = endweek;
			%end;
			%else %if &stop = 1 %then %do;
				stopweek = weektime;
			%end;
		%end;
		%else %if &start = 1 %then %do;
			startweek = weektime;
			stopweek = endweek;
		%end;

		%if &pooled_l = 1 %then %do;
			time = &timevar;
			timesq = time*time;
		%end;


		
	run;

	proc sort data = replic_surv;
	by rand;
	run;

	%if &pooled_l = 0 %then %do;
		
		proc phreg data = replic_surv covm covs(aggregate);
			%if &subset ne . %then %do; by &subset;%end;
			%if &censA = 1 %then %do; 
				where censA = 0;
			%end;
			class sex (ref = "2")  / param = ref;
			%if &class = 1 %then %do;
				class adher (ref = "0")/ param = ref;
			%end;
			%if &class = 2 %then %do;
				class adher (ref = "0") nyha_b (ref = "4")/ param = ref;
			%end;
			%if &ag = 1 %then %do;
				model(startweek, stopweek)*&event(0) = adher &covs /rl; 
			%end;
			%else %do;
				model startweek*&event(0) = adher &covs /rl; 
			%end;
			id &id;
			where startweek < stopweek;
			title "&title_main";

			ods output parameterestimates = param1;
		run;	
	%end;

	%else %if &pooled_l = 1 %then %do;
	
		%if &boots = 1 %then %do;
			
			data params1;
				bsample =.;
			run;

		
			%do bsample = 0 %to &nboot;

				data bootsample;
					merge replic_surv _idsamples (where = (bsample = &bsample));
					by _id;
				run;
				proc sort data = bootsample  sortsize=5G ;
					by b_patid &timevar ;
				run;
				%if &bsample = 0 %then %do;
					
					data bootsample;
						set bootsample;
						numberhits = 1;
					run;
				%end;

				proc logistic data = bootsample descending;
					%if &subset ne . %then %do; by &subset;%end;
					%if &censA = 1 %then %do; 
						where censA = 0;
					%end;
					class sex (ref = "2")  / param = ref;
					%if &class = 1 %then %do;
						class adher (ref = "0")/ param = ref;
					%end;
					%if &class = 2 %then %do;
						class adher (ref = "0") nyha_b (ref = "4")/ param = ref;
					%end;
					model &event = adher time timesq &covs ; 
					title "&title_main";
					freq numberhits;
					
					ods output parameterestimates = param_boots;
				run;

				data param_boots;
				set param_boots;
					bsample = &bsample;
				run;
				
				data params1;
					set params1 param_boots;
					by bsample;
				run;
				proc printto ;
				run;

			%end;		
		%end;
		%else %if &boots= 0 %then %do;

		proc logistic data = replic_surv descending;
			%if &subset ne . %then %do; by &subset;%end;
			%if &censA = 1 %then %do; 
				where censA = 0;
			%end;
			class sex (ref = "2")  / param = ref;
			%if &class = 1 %then %do;
				class adher (ref = "0")/ param = ref;
			%end;
			%if &class = 2 %then %do;
				class adher (ref = "0") nyha_b (ref = "4")/ param = ref;
			%end;
			model &event = adher time timesq &covs ; 
			title "&title_main";
			
			ods output parameterestimates = param1;
		run;

		%end;
	
	%end;


	 %if &boots = 1 %then %do;
		data params2;
		set params1;
		Parameter = Variable;
		if Parameter in ("age0", "hr0", "dbp0", "sbp0") then do; 
			 HazardRatio_new = exp(Estimate*10);
		end;
		else do;
			HazardRatio_new = exp(Estimate);
		end;
		run;
		
		proc sort data = params2; by Parameter;
		run;		

	  	proc univariate data = params2 (where = (bsample >0));
			by Parameter;
			var HazardRatio_new;
			output out = diffpctls (keep = Parameter HR_2_5 HR_97_5)  pctlpre = HR_  pctlpts = 2.5, 97.5;
   		run;
   		data sample0;
			set params2 (where=(bsample = 0));
		    	keep Parameter HazardRatio_new;
			%if &subset ne . %then %do; 
				keep &subset; 
			%end;
		run;

		data pvalues1;
			merge sample0 diffpctls;
			by Parameter;
			HRLowerCL_new = HR_2_5;
			HRUpperCL_new = HR_97_5;	
			if &pooled_l = 1 then model_type = "pooled logistic";			
		run;

		proc print data = pvalues1;
		title 'pvalues1';
		run;

		proc print data = params2;
		title 'params2';
		run;

	%end;


	%else %if &boots = 0 %then %do;
	data pvalues1;
	set param1;

		%if &pooled_l = 1 %then %do;
			Parameter = Variable;
				if Parameter in ("age0", "hr0", "dbp0", "sbp0") then do; 
					 HazardRatio_new = exp(Estimate*10);
					 HRLowerCL_new = exp(Estimate*10 - 1.96*StdErr);
		 			 HRUpperCL_new = exp(Estimate*10 + 1.96*StdErr);
				end;
				else do;
					HazardRatio_new = exp(Estimate);
				 	HRLowerCL_new = exp(Estimate - 1.96*StdErr);
			 		HRUpperCL_new = exp(Estimate + 1.96*StdErr);
				end;
		%end;	

		%else %if &pooled_l = 0 %then %do;
			if StdErrRatio = . then delete;
		
			if Parameter in ("age0", "hr0", "dbp0", "sbp0") then do; 
				 HazardRatio_new = exp(Estimate*10);
				 HRLowerCL_new = exp(Estimate*10 - 1.96*StdErr);
		 		 HRUpperCL_new = exp(Estimate*10 + 1.96*StdErr);
			end;
			else do;
			 	HazardRatio_new = HazardRatio;
				HRLowerCL_new = HRLowerCL;
			 	HRUpperCL_new = HRUpperCL;
			end;
		%end;

		if &pooled_l = 1 then model_type = "pooled logistic";
		else if &pooled_l = 0 then model_type = "cox";

		%if &subset ne . %then %do; 
			keep &subset; 
		%end;
		keep Parameter HazardRatio_new HRLowerCL_new HRUpperCL_new ProbChiSq StdErr  model_type ;
	run;
	%end;

	proc printto print = &outdest;
	run;

	%if &subset ne . %then %do;
		proc sort data = pvalues1;
		by &subset;
		proc print data = pvalues1;
		by &subset; 
		var Parameter Rand HazardRatio_new HRLowerCL_new HRUpperCL_new model_type; 
		format 	&subset rand.;
		title1 &title_main;
		%if &boots = 1 %then %do; title2 "95CI from &nboot bootstrap samples"; %end;
		%else %do; title2 "95% CI from regression"; %end;
		run;
	%end;
	%else %do; 
		proc print data = pvalues1;
		var Parameter HazardRatio_new HRLowerCL_new HRUpperCL_new model_type; 
		title1 &title_main;
		title2 "95% CI from regression";
		run;
	%end;
	proc printto ;
	run;

%end;

%mend nostdz;


/*****************************************************************/
/*Macro stdz: Unadjusted & adjusted analyses with standardization*/
/*****************************************************************/


%macro stdz(outdest = , inset = , id = b_patid, title_main = , nboot= , adjust = , expvar0 = , expvar = , expvar_eof = ,  w_expvar = &expvar,  arm =, eof = ,
	event = , measurevar = , censLTFU = ,  timevar= , visitvar = , recode= 0 , pha = 0, graph = 0,
	itt = 0, ppe = 0, crude = 0, cuminc = 0, ipcw = , iptw = , covs = , covs_t = , weights = stabw1, weight_tests = 0,
	censA = , n_cens = , censvarname =  , protocol = 0, contA = 0, doseform = 0, outc_only = 0, nocens = 0, 
	cens_covs = &covs, cens_covs_t = &covs_t, out_covs = 1, pastA = , trunc_p = 95);

title &title_main;

/*weights indicator*/
%let protocol_ = 0;
%if &iptw = 0 and &ipcw = 0 and &adjust = 1 %then %do;
	%let adjust = 0;
%end;

/*Set up dataset for bootstraps and calculate restricted cubic spline of time*/
proc sort data=&inset out=onesample;
by b_patid &timevar;
run;

%if &itt = 0 and &ppe = 0 %then %do;
data onesample;
set onesample; 
where rand = &arm;
run;
%end;

data onesample ;
  set onesample end = _end_ ;
  by b_patid;
 retain _id ;
  if _n_ = 1 then _id = 0;
  if first.b_patid then do;
  	_id = _id + 1 ;	
  end;
  if _end_ then do ;
     call symput("nids",trim(left(_id)));
  end;

  timesq = &timevar*&timevar;


/*to dichotomize cumulative average of categorical adherence, on natural scale*/
if &recode =1 then do;
	if . < &expvar < 2 then adhvar_t = 1;
	else if  &expvar ge 2 then adhvar_t = 0;

	if . < &expvar0 < 2 then adhvar_0 = 1;
	else if  &expvar0 ge 2 then adhvar_0 = 0;

	%if &cuminc = 1 & &adjust = 1 %then %do;
		if . < &expvar_eof < 2 then adhvar_eof = 1;
		else if  &expvar_eof ge 2 then adhvar_eof = 0;
	%end;		

end;
else if &recode = 2 then do;
	if . < &expvar < 0.8 then adhvar_t = 0;
	else if  &expvar ge 0.8 then adhvar_t = 1;

	if . < &expvar0 < 0.8 then adhvar_0 = 0;
	else if  &expvar0 ge 0.8 then adhvar_0 = 1;

	%if &cuminc = 1 & &adjust = 1 %then %do;
		if . < &expvar_eof < 0.8 then adhvar_eof = 0;
		else if  &expvar_eof ge 0.8 then adhvar_eof = 1;
	%end;	
end;

else if &recode = -1 then do;
	if . < &expvar < 2 then adhvar_t = 1;
	else if  &expvar ge 2 then adhvar_t = 0;
	else if &expvar = . then adhvar_t = 0;

	if . < &expvar0 < 2 then adhvar_0 = 1;
	else if  &expvar0 ge 2 then adhvar_0 = 0;
	else if &expvar0 = . then adhvar_0 = 0;

	%if &cuminc = 1 & &adjust = 1 %then %do;
		if . < &expvar_eof < 2 then adhvar_eof = 1;
		else if  &expvar_eof ge 2 then adhvar_eof = 0;
		else if &expvar_eof = . then adhvar_eof = 0;
	%end;

end;	


else if &recode = -2 then do;
	adhvar_t = 3- &expvar;
	adhvar_0 = 3- &expvar0;

	%if &cuminc = 1 & &adjust = 1 %then %do;
		adhvar_eof = 3 - &expvar_eof;
	%end;
end;	

else if &recode = -3 then do;
	if  &expvar =1 then adhvar_t = 1;
	else if  &expvar = 0 then adhvar_t = 0;
	else if &expvar = . then adhvar_t = 0;

	if  &expvar0 =1 then adhvar_0 = 1;
	else if  &expvar0 = 0 then adhvar_0 = 0;
	else if &expvar0 = . then adhvar_0 = 0;

end;



if &contA = 1 then do;
	if . < &w_expvar < 0.8 then w_adhvar_t = 0;
	else if  &w_expvar ge 0.8 then w_adhvar_t = 1;

	if . < &expvar0 < 0.8 then w_adhvar_0 = 0;
	else if  &expvar0 ge 0.8 then w_adhvar_0 = 1;

	w_adhvar_t1 = lag(w_adhvar_t);

/*Dose response functional forms: 1 = cum(At); 2 = cum(At) + cum(At)^2; 3 = At + cum(At-1); 4 = At + cum(At-1) + cum(At-1_^2)*/
	if &doseform in (1,2) then do;
		adhvar_t = &expvar;
		adhvar_0 = &expvar0;
	 
		if &doseform = 2 then do;
			cumadh_sq = adhvar_t*adhvar_t;
		end;
	end;
	else if &doseform in (3,4) then do;
		adhvar_t = &expvar;
		adhvar_0 = &expvar0;

		cumadh_l1 = lag(adhvar_t);
		if first.b_patid then cumadh_l1 = adhvar_0;
	 
		if &doseform = 4 then do;
			cumadh_l1sq = cumadh_l1*cumadh_l1;
		end;
	end;	
end;



run;

data ids ;
   do bsample = 1 to &nboot;
       do _id = 1 to &nids ;
           output ;
       end;
   end;
run;
proc surveyselect data= ids 
         method = urs
         n= &nids
         seed = 1232  
         out = _idsamples (keep = bsample _id  numberhits  ) 
         outall  noprint  ;       
      strata bsample ;
run;


/*create new dataset with censoring variable*/
proc printto print = &outdest;
run;
proc freq data = onesample;
tables adhvar_0 /list missing;
title 'check 1';
run;
proc printto;
run;

%if &contA = 1 %then %do;
data onesample; 
set onesample; 

censA = 0;

run;

%end;

%else %if &itt = 0 and &ppe = 0 and &contA = 0 and &outc_only ne 1 and &nocens = 0 %then %do;
proc sort data = onesample;
	by b_patid &timevar;

data onesample;
	set onesample;

	by b_patid &timevar;

	if &cuminc = 0 or &adjust = 1 then do;
		retain censA;
	
		if first.b_patid then censA = 0;
		else do;
			if censA = 1 or censA = . then censA = .;
			else if censA = 0 then do;
				if adhvar_t ne adhvar_0 then censA = 1;
				else if adhvar_t = adhvar_0 then censA = 0;
			end;
		end;
		if censA = . then delete;
	end;
run;

%end;
	
%else %if &outc_only = 1 or &nocens = 1 %then %do;
	
	data onesample;
	set onesample;
	
		censA = 0;
	run;
%end;

%else %if &ppe = 1 and &contA = 0 %then %do;
proc sort data = onesample;
	by b_patid &timevar;

	data onesample;
		set onesample;

		by b_patid &timevar;

		retain censA;

		%if &protocol = 1 %then %do;
			if first.b_patid then do;
				if (adhvar_0 = 0 or &censvarname = 1) and everdisc_r ne 1 then censA = 1;
				else censA = 0;
			end;
			else do;
				if censA = 1 or censA = . then censA = .;
				else if censA = 0 then do;
					if (adhvar_t ne 1 or &censvarname = 1) and everdisc_r ne 1 then censA = 1;
					else if adhvar_t = 1 or everdisc_r = 1 then censA = 0;
				end;
			end;
		%end;
		%else %if &protocol = 2 %then %do;
			if first.b_patid then do;
				if (adhvar_0 = 0 or &censvarname = 1) and everdisc_r not in (1, 2) then censA = 1;
				else censA = 0;
			end;
			else do;
				if censA = 1 or censA = . then censA = .;
				else if censA = 0 then do;
					if (adhvar_t ne 1 or &censvarname = 1) and everdisc_r not in (1, 2) then censA = 1;
					else if adhvar_t = 1 or everdisc_r = 1 then censA = 0;
				end;
			end;
		%end;

		%if &adjust = 0 %then %do;
			if censA = . then delete;
		%end;
	run;
%end;

%else %if &itt = 1 %then %do;
data onesample;
set onesample;
	censA = 0;

	adhvar_t = rand;
	adhvar_0 = rand;
run;
%end;


proc printto print = &outdest;
run;
proc freq data = onesample;
tables censA adhvar_t adhvar_0 / missing list;
title 'check 1';
run;
proc printto;
run;


/*create results dataset*/
data means_all;
bsample =.;
&timevar = .;
run;
data hr_all;
bsample =.;
&timevar = .;
run;

%do bsample = 0 %to &nboot;

	title "Bootstrap number &bsample";

	/*set up bootstrap sample*/
	proc sort data = onesample ;
		by _id;
	run;
	data bootsample;
		merge onesample _idsamples (where = (bsample = &bsample));
		by _id;
	run;
	proc sort data = bootsample  sortsize=5G ;
		by b_patid &timevar ;
	run;

	%if &bsample = 0 %then %do;

/*		proc printto print = &outdest;
		run;
*/
		data bootsample;
			set bootsample;
			numberhits = 1;
		run;
	%end;
	
	/*generate IP weights*/
	%if &adjust = 1 %then %do;
		%if &contA = 0 %then %do;
			%weights(datain = bootsample, dataout = trunc, boot = &bsample, weight_tests = &weight_tests,
			ppe=&ppe, iptw = &iptw, ipcw = &ipcw, censvarLTFU = &censLTFU , n_cens = &n_cens, 
			cens_covs = &cens_covs, cens_covs_t = &cens_covs_t, contA = 0, pastA = &pastA, trunc_p = &trunc_p, protocol = &protocol_);	
		%end;
		%else %if &contA = 1 %then %do;
			%weights(datain = bootsample, dataout = trunc, boot = &bsample,  adhvar_t = w_adhvar_t, weight_tests = &weight_tests,
			ppe=&ppe, iptw = &iptw, ipcw = &ipcw, censvarLTFU = &censLTFU , n_cens = &n_cens, 
			cens_covs = &cens_covs, cens_covs_t = &cens_covs_t, contA = 1, pastA = &pastA, trunc_p = &trunc_p, protocol = &protocol_);	

		%end;
		%let datain2 = trunc;
	%end;

	%if &adjust = 0 %then %do;
		%let datain2 = bootsample;
	%end;

	data &datain2;
		set &datain2;

		if &ppe = 1 then do;
			A = rand;
			if censA ne 0 then A = .;
		end;
		else if &ppe = 0 then do;
			A = adhvar_t;
			if &cuminc = 1 and &adjust = 1 then do;
				A = adhvar_eof;
			end;
			if &contA = 1 then do;
				if &doseform in (1, 2) then do;
					A = adhvar_t;
				end;
				else if &doseform in (3, 4) then do;
					A = compyes_t;
				end;
			end;		
		end;
		if &pha = 0 then do;
			Atime = A*&timevar;
			if &doseform = 2 then do;
				cumAsqtime = cumadh_sq*&timevar;
			end;
			else if &doseform = 3 then do;
				cumadhl1time = cumadh_l1*&timevar; 
			end;
			else if &doseform = 4 then do;
				cumadhl1time = cumadh_l1*&timevar; 
				cumAl1sqtime = cumadh_l1sq*&timevar;
			end;
		end;
	run;
	
	%if &pha = 0 %then %do;
		%if &contA = 0 or &doseform = 1 %then %do;
			%let adh_form = A Atime;
		%end;
		%else %if &contA = 1 and &doseform = 2 %then %do;
			%let adh_form = w_adhvar_0  A cumadh_sq Atime cumAsqtime;
		%end;
		%else %if &contA = 1 and &doseform = 3 %then %do;
			%let adh_form = w_adhvar_0  A cumadh_l1 Atime cumadhl1time;	
		%end;
		%else %if &contA = 1 and &doseform = 4 %then %do;
			%let adh_form =  w_adhvar_0 A cumadh_l1 cumadh_l1sq Atime cumadhl1time cumAl1sqtime;
		%end;
	%end;
	%else %if &pha = 1 %then %do;
		%if &contA = 0 or &doseform = 1 %then %do;
			%let adh_form = A;
		%end;
		%else %if &contA = 1 and &doseform = 2 %then %do;
			%let adh_form = w_adhvar_0  A cumadh_sq ;
		%end;
		%else %if &contA = 1 and &doseform = 3 %then %do;
			%let adh_form =  w_adhvar_0 A cumadh_l1;	
		%end;
		%else %if &contA = 1 and &doseform = 4 %then %do;
			%let adh_form = w_adhvar_0  A cumadh_l1 cumadh_l1sq ;
		%end;
	%end;

	%if &cuminc = 1 %then %do;

		proc sort data = &datain2;
		by b_patid &timevar;
		data &datain2;
		set &datain2 (where = (death ne .));
			by b_patid &timevar;
			if last.b_patid then output;
		run;
	%end;

	%if &crude = 0 %then %do;
		%if &cuminc = 1 %then %do;
			/*Run outcome logistic regression model*/
			/*Pr(Y=1|Adherence, Baseline covariates)*/
			proc logistic data = &datain2 descending;
				ods output ParameterEstimates = PE;			       
				model &event =  &adhform  &covs ;
				%if &adjust = 1 %then %do;
					weight &weights;	
				%end;
				title "&title_main , outcome model";
				freq numberhits;
			run;
		%end;
		%else %if &cuminc = 0 %then %do;
			%if &bsample = 0 %then %do;
				proc printto print = &outdest;
				run;
			%end;
			/*Run outcome pooled logistic model*/	
			/*Pr(Yt=1|Adherence, Baseline covariates)*/
			proc logistic data = &datain2 descending;
				where  censA = 0;
				ods output ParameterEstimates = PE;
				%if &out_covs = 0 %then %do;
        				model &event = &adh_form &timevar timesq ;
				%end;
				%else %do;
					model &event = &adh_form &timevar timesq &covs ;
				%end;
				%if &adjust = 1 and (&iptw = 1 or &ipcw = 1) %then %do;
					weight &weights;	
				%end;
				title "&title_main , outcome model";
				freq numberhits;
			run;
			proc freq data = &datain2;
				tables censA /list missing;
				title 'Person-time censored due to non-adherence';
			run;
			proc printto;
			run;
		%end;
	%end;

	%else %if &crude = 1 %then %do;
		%if &cuminc = 1 %then %do;
			/*Run outcome logistic regression model*/
			/*Pr(Y=1|Adherence, Baseline covariates)*/
			proc logistic data = &datain2 descending;
				ods output ParameterEstimates = PE;
			        model &event =  A ;
				freq numberhits;
				title "&title_main , outcome model";
			run;
		%end;
		%else %if &cuminc = 0 %then %do;
			/*Run outcome pooled logistic model*/	
			/*Pr(Yt=1|Adherence, Baseline covariates)*/
			proc logistic data = &datain2 descending;
				where  censA = 0;
				ods output ParameterEstimates = PE;
			  	     model &event = &timevar timesq  &adh_form  ;
				title "&title_main , outcome model";
				freq numberhits;
			run;
		%end;
	%end;
	proc printto ;
	run;

	%if &outc_only ne 0 %then %do;
		data param_boots;
			set PE;
			bsample = &bsample;
		run;
				
		data HR_all;
			set HR_all param_boots;
			by bsample;
		run;

	%end;

	%else %if &outc_only = 0 %then %do;
	/*Using predicted probabilities from appropriate model above, generate standardized dataset and Kaplan-Meier survival estimates*/ 		
		
	proc sql noprint;
		select ESTIMATE FORMAT =16.12 INTO: IBC_ESTIMATE separated by ' ' from pe;
	quit;
	proc sql noprint;
		select variable INTO: model separated by ' ' from PE;
	quit;

	proc means sum noprint data = pe;	
		var df;
		output out = nobs (drop = _type_ _freq_ where=(_stat_ ="N"));
	run;
	proc sql noprint;
		select df into:nvar separated by ' ' from nobs;		
	quit;

	/*create data for three interventions: (1) natural course: &adher = 0, adher = -1, 
	(2) 0% non-adherent: &adher = 1, adher = 0, (3) 100% non-adherent: &adher = 2, adher = 1*/
	/*interpretation note: 0% non-adherent translates to 100% of the time adherent to 
	at least 80% of medication which is the reference level of interest in the analyses*/
 
	%do adher = 0 %to 2;
		%let name_a = rq1;
		%let name = &name_a.&adher;

		data &name (keep = s ci adher &timevar numberhits);		
			set &datain2;
			%if &cuminc = 0 %then %do;
				where &timevar = 0;
			%end;
			array var{&nvar} &model;
			array coef{&nvar} (&ibc_estimate);

			intercept = 1;
			/*numberhits = 1;*/
			adher = &adher - 1;

			/*Expand dataset and calculate predicted survival and risk for natural course*/
			if &adher = 0 then do;	
			    if &cuminc = 0 then do;
				s=1;
				retain A;
				do &timevar = 0 to &eof;
					if &contA = 1 and &doseform = 2 then do;
						cumadh_sq = A*A;
					end;
					else if &contA = 1 and &doseform in (3, 4) then do;
						cumadh_l1 = lag(A);
						if &doseform = 4 then do;
							cumadh_l1sq = cumadh_l1*cumadh_l1sq;
						end;
					end;
					timesq = &timevar*&timevar;
					if &pha = 0 then do;
						Atime = A*&timevar;
						if &contA = 1 and &doseform = 2 then do;
							cumAsqtime = cumadh_sq*&timevar;
						end;
						else if &doseform = 3 then do;
							cumadhl1time = cumadh_l1*&timevar; 
						end;
						else if &doseform = 4 then do;
							cumadhl1time = cumadh_l1*&timevar; 
							cumAl1sqtime = cumadh_l1sq*&timevar;
						end;
					end;
					xbeta = 0;
					do i = 1 to dim(var);
						xbeta = xbeta + coef[i] *var[i];
					end;
        	      			p = 1/(1+exp(-xbeta));
					s = s*(1-p);
					ci = 1-s;
					output;
				end;
			    end;
			    else if &cuminc = 1 then do;
					xbeta = 0;
					do i = 1 to dim(var);
						xbeta = xbeta + coef[i] *var[i];
					end;
        	      			ci = 1/(1+exp(-xbeta));
					&timevar = &eof;
					output;
			    end;
			end;

			/*Expand dataset and calculate predicted survival and risk for interventions*/
			else if &adher in (1, 2) then do;	
			   	A = &adher - 1;
				%if &contA = 1 %then %do;
					 w_adhvar_0  = A;
				%end;
				if &contA = 1 then do;
					if &doseform = 2 then do;
						cumadh_sq = A*A;
					end;
					else if &doseform in (3,4) then do;
						cumadh_l1 = &adher - 1;
						if &doseform = 4 then do;
							cumadh_l1sq = cumadh_l1*cumadh_l1sq;
						end;
					end;
				end;

			    if &cuminc = 0 then do;
				s=1;	
				do &timevar = 0 to &eof;
					if &contA = 1 and &doseform = 2 then do;
						cumadh_sq = A*A;
					end;
					else if &contA = 1 and &doseform in (3, 4) then do;
						cumadh_l1 = lag(A);
						if &doseform = 4 then do;
							cumadh_l1sq = cumadh_l1*cumadh_l1sq;
						end;
					end;
					timesq = &timevar*&timevar;
					if &pha = 0 then do;
						Atime = A*&timevar;
						if &contA = 1 and &doseform = 2 then do;
							cumAsqtime = cumadh_sq*&timevar;
						end;
						else if &doseform = 3 then do;
							cumadhl1time = cumadh_l1*&timevar; 
						end;
						else if &doseform = 4 then do;
							cumadhl1time = cumadh_l1*&timevar; 
							cumAl1sqtime = cumadh_l1sq*&timevar;
						end;
					end;
					xbeta = 0;
					do i = 1 to dim(var);
						xbeta = xbeta + coef[i] *var[i];
					end;	
        		      		p = 1/(1+exp(-xbeta));
					s = s*(1-p);
					ci = 1-s;
					output;
				end;
			    end;
    			    else if &cuminc = 1 then do;
					xbeta = 0;
					do i = 1 to dim(var);
						xbeta = xbeta + coef[i] *var[i];
					end;	
        		      		ci = 1/(1+exp(-xbeta));
					&timevar = &eof;
					output;
			    end;
			end;	
		run;
	%end;

	/*combine all three interventions*/
	data rq1;
		set rq10 rq11 rq12;
		by adher;
	run;

	/*calculate mean cumulative incidence for each visit and intervention*/
	proc means data = rq1 mean noprint;
		class &timevar adher;
		types &timevar*adher;
		var ci;
		freq numberhits;
		output out = mean_1 (drop =_type_ _freq_) mean(ci) = ;
	run;
	data mean_1;
		set mean_1;
		label ci = "Cumulative incidence" ;
		bsample = &bsample;
	run;

	/*combine across bsamples*/
	data means_all;
		set means_all mean_1;
		by bsample &timevar;
		if bsample = . then delete;
	run;

	%if &cuminc = 0 %then %do;
		data hr;
	 	set pe;	
			%if &contA = 1 %then %do;
				where Variable in ('A', 'cumadh_l1', 'cumadh_l1sq');
			%end;
			%else %do;
				where Variable = 'A';
			%end;
			bsample = &bsample;
			
		keep bsample variable estimate;
		run;
	
	data hr_all;
		set hr_all hr;
		by bsample; 
		if bsample = . then delete;
	run;

	%end;
	proc datasets library = work nolist;
		delete censadh_num0  censadh_dnom0 temp pctl rq1 pe mean_1;
	run;
	proc printto ;
	run;

  %end;

%end;

%if &outc_only = 0 %then %do;

proc printto print = &outdest;
run;
title "Summary";
proc sort data=means_all;
by bsample &timevar;
run;

/*Calculate standardized 5-year risk difference and 95% confidence interval*/
/*Comparing Always adhere to less than 80% to Always adhere to at least 80% of medication (adher: 1 vs 0)*/
proc transpose data=means_all out = temp prefix = Risk_;
var ci;
id adher;
by bsample &timevar;
run;

data temp;
set temp;
  	rd = Risk_1 - Risk_0;
	surv_0 = 1-risk_0;
	surv_1 = 1-risk_1;
	%if &cuminc = 0 %then %do;
		ratio = log(surv_1)/log(surv_0);
	%end;
	%if  &cuminc = 1 %then %do;
		ratio = (risk_1/(1-risk_1))/(risk_0/(1-risk_0));
	%end;
run;

%if &cuminc = 0 %then %do;
	proc sort data = temp;
	by bsample;
	run;
	proc means data = temp mean noprint ;
	by bsample;
	var ratio;
	output out= temp2 (keep = mean_ratio bsample) mean = mean_ratio;
	run;

	data temp;
	merge temp temp2;
	by bsample;
		ratio_t = ratio;
		ratio = mean_ratio;
	run;
%end;

%if &cuminc = 1 %then %do;
	data temp;
	set temp;
	ratio_t = .;
	run;
%end;

proc sort data = temp;
by &timevar;
   proc univariate data = temp (where = (bsample >0)) noprint;
 	by &timevar;
	var rd risk_0 risk_1 surv_0 surv_1 ratio_t ratio ;
	output out = diffpctls  pctlpre = rd_ p0_ p1_ s0_ s1_ ratio_t_ ratio_  pctlpts = 2.5, 97.5;
   run;
   data sample0;
	set temp (where=(bsample = 0));
     	keep  rd  risk_0 risk_1 surv_0 surv_1 ratio_t ratio &timevar;
   run;


proc sort data= diffpctls;  by &timevar;
run;

data final;
	merge sample0 diffpctls; 
	by &timevar;
	
	%if &cuminc = 0 %then %do;
		avgHR = ratio;
		avgHR_2_5 = ratio_2_5;	
		avgHR_97_5 = ratio_97_5;

		cHR = ratio_t;
		cHR_2_5 = ratio_t_2_5;	
		cHR_97_5 = ratio_t_97_5;

	%end;

	%if &cuminc = 1 %then %do;
		OR = ratio;
		OR_2_5 = ratio_2_5;
		OR_97_5 = ratio_97_5;
	%end;
run;


/*output final results to .rtf file*/

proc printto print = &outdest;
run;

%if &pha = 1 and &cuminc = 0 %then %do;
	proc print data = final label noobs;
	where &timevar in (0, 26, 52, 78, 104, 130, 156, 182, 205);
		var &timevar risk_0 risk_1 p0_2_5 p0_97_5 p1_2_5 p1_97_5 ;
	title1 "&title_main";
	title2 "Risk and risk difference";
	title3 "95% Confidence Intervals using &nboot samples" ;
	%if &arm = 0 %then %do; title4  'Placebo';%end; %else %if &arm = 1 %then %do; title4 'Candesartan'; %end; %else %if &arm = . %then %do; title4 'Overall';%end;
	run;

	proc print data = final label noobs;
	where &timevar in (0, 26, 52, 78, 104, 130, 156, 182, 205);
		var &timevar surv_0 surv_1 s0_2_5 s0_97_5 s1_2_5 s1_97_5 /*avgHR avgHR_2_5 avgHR_97_5 */ cHR cHR_2_5 cHR_97_5; 
	title2 "Survival and HR";
	title1 "&title_main";
	title3 "95% Confidence Intervals using &nboot samples" ;
	%if &arm = 0 %then %do; title4  'Placebo';%end; %else %if &arm = 1 %then %do; title4 'Candesartan'; %end; %else %if &arm = . %then %do; title4 'Overall';%end;
	run;
%end;
%if &cuminc = 1 %then %do;
	proc print data = final label noobs;
	where &timevar in (0, 26, 52, 78, 104, 130, 156, 182, 205);
	var risk_0 risk_1 p0_2_5 p0_97_5 p1_2_5 p1_97_5 rd  rd_2_5 rd_97_5;

	proc print data = final label noobs;
	where &timevar in (0, 26, 52, 78, 104, 130, 156, 182, 205);
	var OR OR_2_5 OR_97_5; 
	title2 "Odds Ratio";
	title1 "&title_main";
	title3 "95% Confidence Intervals using &nboot samples" ;
	%if &arm = 0 %then %do; title4  'Placebo';%end; %else %if &arm = 1 %then %do; title4 'Candesartan'; %end; %else %if &arm = . %then %do; title4 'Overall';%end;
	run;
%end;
%if &pha = 0 %then %do;
	proc print data = final label noobs;
	where &timevar in (0, 26, 52, 78, 104, 130, 156, 182, 205);
		var &timevar risk_0 risk_1 p0_2_5 p0_97_5 p1_2_5 p1_97_5 rd  rd_2_5 rd_97_5; 
	title1 "&title_main";
	title2 "Risk and risk difference";
	title3 "95% Confidence Intervals using &nboot samples" ;
	%if &arm = 0 %then %do; title4  'Placebo';%end; %else %if &arm = 1 %then %do; title4 'Candesartan'; %end; %else %if &arm = . %then %do; title4 'Overall';%end;
	run;
	proc print data = final label noobs;
	where &timevar in (0, 26, 52, 78, 104, 130, 156, 182, 205);
	var &timevar surv_0 surv_1 s0_2_5 s0_97_5 s1_2_5 s1_97_5  ; 
	title2 "Survival ";
	run;

	%if &graph = 1 %then %do;
		proc print data = final label noobs;
			var &timevar risk_0 risk_1 p0_2_5 p0_97_5 p1_2_5 p1_97_5 rd  rd_2_5 rd_97_5; 
		title2 "Cumulative incidence curves";
		title1 "&title_main";
		title3 "95% Confidence Intervals using &nboot samples" ;
		%if &arm = 0 %then %do; title4  'Placebo';%end; %else %if &arm = 1 %then %do; title4 'Candesartan'; %end; %else %if &arm = . %then %do; title4 'Overall';%end;
		run;

		proc print data = final label noobs;
			var &timevar surv_0 surv_1 s0_2_5 s0_97_5 s1_2_5 s1_97_5 ; 
		title2 "Survival curves";
		title1 "&title_main";
		title3 "95% Confidence Intervals using &nboot samples" ;
		%if &arm = 0 %then %do; title4  'Placebo';%end; %else %if &arm = 1 %then %do; title4 'Candesartan'; %end; %else %if &arm = . %then %do; title4 'Overall';%end;
		run;
	%end;
%end;

%end; 

%else %if &outc_only ne 0 %then %do;

	data params2;
		set HR_all;

		if Variable in ("age0", "hr0", "dbp0", "sbp0") then do; 
			 HazardRatio_new = exp(Estimate*10);
		end;
		else do;
			HazardRatio_new = exp(Estimate);
		end;
		run;
		
		proc sort data = params2; by Variable;
		run;		

	  	proc univariate data = params2 (where = (bsample >0));
			by Variable;
			var HazardRatio_new;
			output out = diffpctls (keep = Variable HR_2_5 HR_97_5)  pctlpre = HR_  pctlpts = 2.5, 97.5;
   		run;
   		data sample0;
			set params2 (where=(bsample = 0));
		    	keep Variable HazardRatio_new;
		run;

		data pvalues1;
			merge sample0 diffpctls;
			by Variable;
			HRLowerCL_new = HR_2_5;
			HRUpperCL_new = HR_97_5;	
			model_type = "pooled logistic";			
		run;

		proc print data = HR_all;
		title 'HR_all';
		run;

		proc print data = params2;
		title 'params2';
		run;

	proc printto print = &outdest;
	run;

	title "Adherence modeled as &expvar = &title_main";

	proc print data = pvalues1;
		var Variable HazardRatio_new HRLowerCL_new HRUpperCL_new model_type; 
		title1 "&title_main";
		title2 "95CI from &nboot bootstrap samples"; 
	run;
	proc printto ;
	run;


%end;

proc printto;
run;
%let timenow2=%sysfunc(time(), time.);
%let datenow2=%sysfunc(date(), date9.);
%put Part C is complete;
%put End time is &datenow2 / &timenow2 ;
%put Program is complete;
%put ; 	
%mend stdz;


/*********************************************************/
/*Macro weights: Inverse probability weights calculation */
/*********************************************************/

%macro weights(datain = , dataout = , boot = , ppe =, iptw = , ipcw = , censvarLTFU = , n_cens = );

	proc printto;
	run;
		
	%if &ppe = 1 %then %do; %let loops = 2; %end;
	%else %if &ppe = 0 %then %do; %let loops = 1; %end;

	%do loop = 0 %to (&loops-1);

		data weights;
			set &datain;
	
			%if &ppe = 1 %then %do; 
				where rand = &loop;
			%end;


			%if &weights = 'stabw1' %then %do;
				array covs(*)  adhvar_t adhvar_0  &covs &covs_t;
			%end;
			%else %do;
				array covs(*)   adhvar_t adhvar_0  &covs;		 
			%end;
			do i = 1 to dim(covs);
				if covs(i) = . then delete;
			end;


		run;

		%if &ipcw = 1 %then %do;
			proc sort data = weights;
				by b_patid newvisit;
			data visits;
			set  weights (where =(&visitvar = 1));
			by b_patid newvisit;

				%if &n_cens ge 2 %then %do;
					measure_l1 = lag(&measurevar);
					if newvisit le 1 then do;
						measure_l1 = 1;
					end;
				%end;
				%else %do;
					measure_l1 = 0; /*set to zero so that lag 2 data not a restriction in weights*/
				%end;
				%if &n_cens ge 3 %then %do;
					measure_l2 = lag2(&measurevar);
					if newvisit le 2 then do;
						measure_l2 = 1;
					end;
				%end;
				%else %do;
					measure_l2 = 0; /*set to zero so that lag 2 data not a restriction in weights*/
				%end;
				%if &n_cens ge 4 %then %do;
					measure_l3 = lag3(&measurevar);
					if newvisit le 3 then do;
						measure_l3 = 1;
					end;
				%end;
				%else %do;
					measure_l3 = 0; /*set to zero so that lag 2 data not a restriction in weights*/
				%end;

				if first.b_patid then do;
					measure_l1 = 1;
					measure_l2 = 1;
					measure_l3 = 1;
				end;
			
			run;
			data  weights;
			merge  weights visits;
			by b_patid newvisit;
			run;

		%end;
		%if &iptw = 1 %then %do;

		  	/*IPW for adherence at time t - point estimate*/

			/*model for censoring due to adherence: censA = 0 vs 1; only eligible for censoring when week is true visit week*/
			/*Numerator: Pr(At = 1|A_0, Visit & Measured, Baseline covariates)*/

			proc logistic data = weights (where =(&visitvar = 1 and  &measurevar = 1)) descending;
				model adhvar_t= adhvar_0 &timevar timesq &covs ;
				freq numberhits;	
				output out = adhcens_num&loop (keep=b_patid weektime adhnum ) p = adhnum;
			run;

			/*Denominator:  Pr(At = 1|A_0, Visit & Measured, Baseline covariates, Time-varying covariates)*/
			proc logistic data = weights (where =(&visitvar = 1 and  &measurevar = 1)) descending;
				model adhvar_t= adhvar_0  &timevar timesq &covs &covs_t;
				freq numberhits;
				output out = adhcens_dnom&loop (keep=b_patid weektime adhdnom ) p = adhdnom;
			run;
		%end;
		%if &ipcw = 1 %then %do;
 
			/*model for loss to follow-up censoring at time t*/
			/*Numerator: Pr(CensLTFU =0|A_0, Baseline covariates, no adherence at t-1 or t-2)*/
			/*comp_m_t1 = 0 if no measurement at t-1 most recent scheduled visit; comp_m_t2 = 0 if no measurement at t-2 most recent scheduled visit*/

			proc logistic data= weights (where =(&visitvar = 1 and measure_l1 = 0 and measure_l2 = 0));
				model &censvarLTFU = &timevar timesq adhvar_0 &covs ;
				freq numberhits;
				output out = censLTFU_num&loop (keep=b_patid weektime censnum ) p =censnum;
			run;

			/*Denominator: Pr(CensLTFU =0|A_0, Baseline covariates, Time-varying covariates, no adherence at t-1 or t-2)*/

			proc logistic data= weights (where =(&visitvar = 1 and measure_l1 = 0 and measure_l2 = 0));
				model  &censvarLTFU = adhvar_0 &timevar timesq	&covs &covs_t ;
					freq numberhits;
					output out = censLTFU_dnom&loop (keep = b_patid weektime censdnom )  p = censdnom;
				run;
		%end;
		%if &iptw = 1 %then %do;
			proc sort data=adhcens_num&loop;
				by b_patid weektime;
			proc sort data=adhcens_dnom&loop;
				by b_patid weektime;
			run;
		%end;
		%if &ipcw = 1 %then %do;
			proc sort data=censLTFU_num&loop;
				by b_patid weektime;
			proc sort data=censLTFU_dnom&loop;
				by b_patid weektime;
			run;
		%end;	

	%end;


		%if &ppe = 0 %then %do;
			%if &iptw = 1 and &ipcw = 1 %then %do;
				%let datasets = adhcens_num0 adhcens_dnom0 censLTFU_num0 censLTFU_dnom0 ;
			%end;
			%else %if &iptw = 1 and &ipcw = 0 %then %do;
				%let datasets = adhcens_num0 adhcens_dnom0 ;
			%end;
			%else %if &iptw = 0 and &ipcw = 1 %then %do;
				%let datasets = censLTFU_num0 censLTFU_dnom0 ;
			%end;
		%end;
		%else %if &ppe = 1 %then %do;			
			%if &iptw = 1 and &ipcw = 1 %then %do;
				%let datasets = adhcens_num0 adhcens_dnom0 censLTFU_num0 censLTFU_dnom0 
				      		adhcens_num1 adhcens_dnom1 censLTFU_num1 censLTFU_dnom1  ;
			%end;
			%else %if &iptw = 1 and &ipcw = 0 %then %do;
				%let datasets = adhcens_num0 adhcens_dnom0 adhcens_num1 adhcens_dnom1; 
			%end;
			%else %if &iptw = 0 and &ipcw = 1 %then %do;
				%let datasets = censLTFU_num0 censLTFU_dnom0 censLTFU_num1 censLTFU_dnom1 ;
			%end;
		%end;		

	proc sort data= &datain;
		by b_patid weektime;
	data main_w6;
		merge &datain &datasets ;
	        by b_patid  weektime;

		/* variables ending with _0 refer to the numerator of the weights
		   Variables ending with _w refer to the denominator of the weights */

		        if first.b_patid then do; 
				k1_0=1;
				k1_w=1; 
				m1_0=1;
				m1_w=1;
			end;
			retain k1_0 k1_w m1_0 m1_w;
		if &iptw = 1 then do;
			if censA = 0 then do; /*not censored due to adherence change*/
				if adhvar_t = 1 then do; /*adherent at t & since baseline*/
					if adhnum  = . then adhnum  = 1;	
					if adhdnom = . then adhdnom = 1;
					m1_0=m1_0*adhnum ;
					m1_w=m1_w*adhdnom ;
				end;
				else if adhvar_t = 0 then do; /*non-adherent at t*/
					if adhnum = . then adhnum = 0;	
					if adhdnom = . then adhdnom = 0;
					m1_0=m1_0*(1-adhnum);
					m1_w=m1_w*(1-adhdnom);
				end;
			end;
			else if censA ne 0 then do; /*censored due to adherence*/
				adhnum = 1;	
				adhdnom = 1;
				m1_0=m1_0*adhnum;
				m1_w=m1_w*adhdnom;
			end;
		end;
		else if &iptw = 0 then do; /*do not use treatment weights*/
			adhnum = 1;	
			adhdnom = 1;
			m1_0=1;
			m1_w=1;
		end;
		if &ipcw = 1 then do; /*use loss to follow-up weights*/
			/*prob. not censored due to loss to follow-up*/
			if censnum = . then censnum =1;
			if censdnom = . then censdnom =1;
			k1_0=k1_0*(censnum);
		        k1_w=k1_w*(censdnom);
		end;
		else if &ipcw = 0 then do; /*do not use loss to follow-up weights*/
			censnum =1;
			censdnom =1;
		        k1_0=1;
		       	k1_w=1;
		end;

		stabw_A=(m1_0)/(m1_w);
		nstabw_A=1/(m1_w);
		
		if k1_w ne 0 then do;		
			stabw_LTFU=(k1_0)/(k1_w);
			nstabw_LTFU=1/(k1_w);
		end;
		else do;
			stabw_LTFU=.;
			nstabw_LTFU=.;
		end;

		stabw =stabw_A*stabw_LTFU;
		nstabw =nstabw_A*nstabw_LTFU;
	run;
	%if &boot = 0 %then %do;
		/*proc printto print = &outdest;
		run;*/
		proc sort data = main_w6; by rand;run;
		proc means data=main_w6  n mean std min max p95 p99 nmiss;
			by rand;
			var nstabw stabw nstabw_a stabw_a nstabw_ltfu stabw_ltfu;
			title 'weights, all';
		run;
		/*proc printto print = &outdest;
		run;*/
		proc freq data=main_w6 nlevels;
			where stabw ne .; 
			by rand;
			tables b_patid /noprint;
			title 'weights, all';
		run;
	%end;
		
	/*for truncation*/
	proc means data=main_w6  n mean std min max p95 p99 nmiss noprint;
		var stabw;
		title 'stabilized weights, end of follow-up';
		output out=pctl (keep =  p99) p99 = p99;
	run;
	proc means data=main_w6   p99 noprint;

		var nstabw;
		title 'stabilized weights, end of follow-up';
		output out=pctl_n (keep = p99) p99 = p99;
	run;	

	data temp;
		set pctl;
		call symput ('cutoff', p99);
		
	run;
	data temp_n;
		set pctl_n;
		call symput ('cutoff_n', p99);
	run;

	data &dataout; 
		set main_w6;
		stabw1 = stabw;
		nstabw1 = nstabw;
		if stabw >  %sysevalf(&cutoff)  then do;
			stabw1 = %sysevalf(&cutoff);
		end;
		if nstabw >  %sysevalf(&cutoff_n)  then do;
			nstabw1 = %sysevalf(&cutoff_n);
		end;
	run;	

	%if &boot = 0 %then %do;
		proc printto print = &outdest;
		run;
		proc means data=&dataout  n mean std min max p95 p99 nmiss;
			by rand;
			var stabw stabw1 nstabw nstabw1;
			title 'stabilized weights, end of follow-up';
		run;
		proc printto;
		run;
	%end;

%mend weights;


/*********************************************************************************/
/*Macro dose_r: Unadjusted and adjusted analyses with dose-response of adherence */
/*********************************************************************************/


%macro dose_r(outdest = , inset = , title_main = , nboot= , adjust = , expvar0 = , expvar = ,expvar_eof = , eof = ,  w_expvar = , 
	event = , measurevar = , censLTFU = ,  timevar= , visitvar = , recode= 2 , pha = 0, graph = 0, weighted = 1,weight_tests = 0,
	ipcw = , iptw = , covs = , covs_t = , weights = stabw1, censA = , n_cens = , censvarname = , protocol =  , doseform = 0, 
	cens_covs = &covs, cens_covs_t = &covs_t, out_covs = 1, pastA = ,  trunc_p = 95, intA = 1);

title &title_main;

/*weights indicator*/
%if &protocol ne 0 %then %do;
	%let protocol_ = &protocol;
%end;
%else %if &protocol = 0 %then %do;
	%let protocol_ = 5;
%end;


/*Set up dataset for bootstraps and calculate restricted cubic spline of time*/
proc sort data=&inset out=onesample;
by b_patid weektime;
run;

data onesample ;
  set onesample end = _end_ ;
  by b_patid;
 retain _id ;
  if _n_ = 1 then _id = 0;
  if first.b_patid then do;
  	_id = _id + 1 ;	
  end;
  if _end_ then do ;
     call symput("nids",trim(left(_id)));
  end;

  timesq = weektime*weektime;
  timecb = timesq*weektime;

  censA = 0;

  /*to dichotomize cumulative average of categorical adherence, for weights*/
	if . < &w_expvar  < 0.8 then w_adhvar_t = 0;
	else if  &w_expvar ge 0.8 then w_adhvar_t = 1;

	if . < &expvar0 < 0.8 then w_adhvar_0 = 0;
	else if  &expvar0 ge 0.8 then w_adhvar_0 = 1;

	w_adhvar_t1 = lag(w_adhvar_t);

	/*Dose response functional forms: 1 = cum(At); 2 = cum(At) + cum(At)^2; 3 = At + cum(At-1); 4 = At + cum(At-1) + cum(At-1_^2)*/
	if &doseform in (1,2) then do;
		adhvar_t = &expvar;
		adhvar_0 = &expvar0;
	 
		if &doseform = 2 then do;
			cumadh_sq = adhvar_t*adhvar_t;
		end;
	end;
	else if &doseform in (3,4) then do;
		adhvar_t = &expvar;
		adhvar_0 = &expvar0;

		cumadh_l1 = lag(adhvar_t);
		if first.b_patid then cumadh_l1 = adhvar_0;
	 
		if &doseform = 4 then do;
			cumadh_l1sq = cumadh_l1*cumadh_l1;
		end;
	end;	
run;

/*
proc print data = onesample (obs = 600);
var b_patid weektime visit_t newvisit adhvar_0 adhvar_t &expvar everdisc_r recode_flag_t cumadh_l1 cumadh_l1sq;
run;
*/

data ids ;
   do bsample = 1 to &nboot;
       do _id = 1 to &nids ;
           output ;
       end;
   end;
run;
proc surveyselect data= ids 
         method = urs
         n= &nids
         seed = 1232  
         out = _idsamples (keep = bsample _id  numberhits  ) 
         outall  noprint  ;       
      strata bsample ;
run;

/*create results dataset*/
data means_all;
bsample =.;
weektime = .;
run;
data hr_all;
bsample =.;
weektime = .;
run;

%do bsample = 0 %to &nboot;

	title "Bootstrap number &bsample";

	/*set up bootstrap sample*/
	proc sort data = onesample ;
		by _id;
	run;
	data bootsample;
		merge onesample _idsamples (where = (bsample = &bsample));
		by _id;
	run;
	proc sort data = bootsample  sortsize=5G ;
		by b_patid weektime ;
	run;

	%if &bsample = 0 %then %do;

		proc printto print = &outdest;
		run;

		data bootsample;
			set bootsample;
			numberhits = 1;
		run;
	%end;
	
	%if &weighted = 1 %then %do;
		/*generate IP weights*/
		%weights(datain = bootsample, dataout = trunc, boot = &bsample, adhvar_t = w_adhvar_t,weight_tests = &weight_tests, trunc_p = &trunc_p,
		ppe=1, iptw = &iptw, ipcw = &ipcw, censvarLTFU = &censLTFU , n_cens = &n_cens, cens_covs = &cens_covs, cens_covs_t = &cens_covs_t, contA = 1, pastA = &pastA, 
		dose_r = 1, protocol = &protocol_);	
		%let datain2 = trunc;
	%end;
	%else %do;
		%let datain2 = bootsample;
	%end;
	data &datain2;
		set &datain2;
		A = adhvar_t;
	
		if &pha = 0 then do;
			Atime = A*&timevar;
			if &doseform = 2 then do;
				cumAsqtime = cumadh_sq*&timevar;
			end;
			else if &doseform = 3 then do;
				cumadhl1time = cumadh_l1*&timevar; 
			end;
			else if &doseform = 4 then do;
				cumadhl1time = cumadh_l1*&timevar; 
				cumAl1sqtime = cumadh_l1sq*&timevar;
			end;
		end;
	run;

	/*Run outcome pooled logistic model*/	
	proc sort data = &datain2;
	by rand;
	run;
	
			
	%if &pha = 0 %then %do;
		%if &doseform = 1%then %do;
			%let adh_form = A Atime;
		%end;
		%else %if &doseform = 2 %then %do;
			%let adh_form = A cumadh_sq Atime cumAsqtime;
		%end;
		%else %if &doseform = 3 %then %do;
			%let adh_form = A cumadh_l1 Atime cumadhl1time;	
		%end;
		%else %if &doseform = 4 %then %do;
			%let adh_form = A cumadh_l1 cumadh_l1sq Atime cumadhl1time cumAl1sqtime;
		%end;
	%end;
	%else %if &pha = 1 %then %do;
		%if &doseform = 1%then %do;
			%let adh_form = A;
		%end;
		%else %if &doseform = 2 %then %do;
			%let adh_form = A cumadh_sq ;
		%end;
		%else %if &doseform = 3 %then %do;
			%let adh_form = A cumadh_l1;	
		%end;
		%else %if &doseform = 4 %then %do;
			%let adh_form = A cumadh_l1 cumadh_l1sq ;
		%end;
	%end;


	%do loop = 0 %to 1;
		%if &bsample = 0 %then %do;
			proc printto print = &outdest;
			run;
		%end;

		
		/*Pr(Yt=1|Adherence, Baseline covariates)*/
		proc logistic data = &datain2 descending;
			ods output ParameterEstimates = PE;
			where rand = &loop;
			%if &out_covs = 0 %then %do;
	        		model &event = &timevar timesq w_adhvar_0  &adh_form ;
			%end;
			%else %do;
	        		model &event = &timevar timesq w_adhvar_0 &adh_form &covs ;
			%end;
			%if &weighted = 1 %then %do;
				weight &weights;	
			%end;
			freq numberhits;
		run;

		proc printto;
		run;
		/*Using predicted probabilities from appropriate model above, generate standardized dataset and Kaplan-Meier survival estimates*/ 		
		
		proc sql noprint;
			select ESTIMATE FORMAT =16.12 INTO: IBC_ESTIMATE separated by ' ' from pe;
		quit;
		proc sql noprint;
			select variable INTO: model separated by ' ' from PE;
		quit;

		proc means sum noprint data = pe;	
			var df;
			output out = nobs (drop = _type_ _freq_ where=(_stat_ ="N"));
		run;
		proc sql noprint;
			select df into:nvar separated by ' ' from nobs;		
		quit;

		/*create data for 100% adherence */
		%let name_a = rq1;
		%let name = &name_a.&loop;

		data &name (keep = s ci loop rand weektime numberhits);		
			set &datain2;
			where weektime = 0;
			array var{&nvar} &model;
			array coef{&nvar} (&ibc_estimate);

			intercept = 1;
			numberhits = 1;
			loop = &loop;
	
			/*Expand dataset and calculate predicted survival and risk for natural course*/
			w_adhvar_0 = 1;
			A = &intA;

			if  &doseform = 2 then do;
				cumadh_sq = A*A;
			end;
			else if &doseform in (3, 4) then do;
				cumadh_l1 = lag(A);
				if &doseform = 4 then do;
					cumadh_l1sq = cumadh_l1*cumadh_l1;
				end;
			end;
			
			s=1;	
			do weektime = 0 to &eof;
				timesq = weektime*weektime;
				if &pha = 0 then do;
					Atime = A*&timevar;
					if &doseform = 2 then do;
						cumAsqtime = cumadh_sq*&timevar;
					end;
					else if &doseform = 3 then do;
						cumadhl1time = cumadh_l1*&timevar; 
					end;
					else if &doseform = 4 then do;
						cumadhl1time = cumadh_l1*&timevar; 
						cumAl1sqtime = cumadh_l1sq*&timevar;
					end;
				end;

				xbeta = 0;
				do i = 1 to dim(var);
					xbeta = xbeta + coef[i] *var[i];
				end;	
   		      		p = 1/(1+exp(-xbeta));
				s = s*(1-p);
				ci = 1-s;
				output;
			end;
		
		run;	
		
		/*calculate mean cumulative incidence for each visit and intervention*/
		proc means data = &name mean noprint;
			class weektime ;
			types weektime;
			var ci;
			freq numberhits;
			output out = mean_1 (drop =_type_ _freq_) mean(ci) = ;
		run;
		data mean_1;
			set mean_1;
			label ci = "Cumulative incidence" ;
			bsample = &bsample;
			rand_new = &loop;
		run;

		/*combine across bsamples*/
		data means_all;
			set means_all mean_1;
			by bsample weektime;
			if bsample = . then delete;
		run;
	%end;		

%end;

proc printto print = &outdest;
run;
title "Summary";
proc sort data=means_all;
by bsample weektime rand_new ;
run;

/*Calculate standardized 5-year risk difference and 95% confidence interval*/
/*Comparing Always adhere to less than 80% to Always adhere to at least 80% of medication (adher: 1 vs 0)*/
proc transpose data=means_all out = temp prefix = Risk_;
var ci;
id rand_new;
by bsample weektime;
run;

data temp;
set temp;

 	rd = Risk_1 - Risk_0;
	surv_0 = 1-risk_0;
	surv_1 = 1-risk_1;
	
	HR = log(surv_1)/log(surv_0);
	
run;


proc sort data = temp;
by bsample;
run;
proc means data = temp mean noprint;
by bsample;
var HR;
output out= temp2 (keep = mean_HR bsample) mean = mean_HR;
run;

data temp;
merge temp temp2;
by bsample;
	cHR =HR;
	HR = mean_HR;
run;


proc sort data = temp;
by weektime;
   proc univariate data = temp (where = (bsample >0)) noprint;
 	by weektime;
	var rd risk_0 risk_1 surv_0 surv_1 HR cHR ;
	output out = diffpctls  pctlpre = rd_ p0_ p1_ s0_ s1_ HR_ cHR_ pctlpts = 2.5, 97.5;
   run;
   data sample0;
	set temp (where=(bsample = 0));
     	keep  rd  risk_0 risk_1 surv_0 surv_1 HR cHR weektime;
   run;


proc sort data= diffpctls;  by weektime;
run;

data final;
	merge sample0 diffpctls; 
	by weektime;
	
run;

/*output final results to .rtf file*/

proc printto print = &outdest;
run;


%if &pha = 1 %then %do;
	proc print data = final label noobs;
	where weektime in (0, 26, 52, 78, 104, 130, 156, 182, 205);
		var weektime risk_0 risk_1 p0_2_5 p0_97_5 p1_2_5 p1_97_5 ;
	title1 "&title_main";
	title2 "Risk and risk difference";
	title3 "95% Confidence Intervals using &nboot samples" ;
	title4 'Overall';
	run;

	proc print data = final label noobs;
	where weektime in (0, 26, 52, 78, 104, 130, 156, 182, 205);
		var weektime surv_0 surv_1 s0_2_5 s0_97_5 s1_2_5 s1_97_5 /*HR HR_2_5 HR_97_5*/ cHR cHR_2_5 cHR_97_5; 
	title2 "Survival and HR";
	title1 "&title_main";
	title3 "95% Confidence Intervals using &nboot samples" ;
	title4 'Overall';
	run;
%end;


%if &pha = 0 %then %do;
	proc print data = final label noobs;
	where weektime in (0, 26, 52, 78, 104, 130, 156, 182, 205);
		var weektime risk_0 risk_1 p0_2_5 p0_97_5 p1_2_5 p1_97_5 rd  rd_2_5 rd_97_5; 
	title1 "&title_main";
	title2 "Risk and risk difference";
	title3 "95% Confidence Intervals using &nboot samples" ;
	title4 'Overall';
	run;
	proc print data = final label noobs;
	where weektime in (0, 26, 52, 78, 104, 130, 156, 182, 205);
	var weektime surv_0 surv_1 s0_2_5 s0_97_5 s1_2_5 s1_97_5  ; 
	title2 "Survival ";
	run;

	%if &graph = 1 %then %do;
		proc print data = final label noobs;
			var weektime risk_0 risk_1 p0_2_5 p0_97_5 p1_2_5 p1_97_5 rd  rd_2_5 rd_97_5; 
		title2 "Cumulative incidence curves";
		title1 "&title_main";
		title3 "95% Confidence Intervals using &nboot samples" ;
		title4 'Overall';
		run;

		proc print data = final label noobs;
			var weektime surv_0 surv_1 s0_2_5 s0_97_5 s1_2_5 s1_97_5 ; 
		title2 "Survival curves";
		title1 "&title_main";
		title3 "95% Confidence Intervals using &nboot samples" ;
		title4 'Overall';
		run;
	%end;
%end;



proc printto;
run;
%let timenow2=%sysfunc(time(), time.);
%let datenow2=%sysfunc(date(), date9.);
%put Part C is complete;
%put End time is &datenow2 / &timenow2 ;
%put Program is complete;
%put ; 	


%mend dose_r;


/**************************************************/
/*Macro Table 1: Intention-to-treat***************/
/*************************************************/
%macro Table1(nboot = , eof = , row = );
	
	%if &row le 2 %then %do;
		/*Cox model*/
		data itt;
		set master.master;
		keep b_patid rand study_ t2dthc death dgbb_b dgdiu_b dgacei_b dglild_b 
			hxaf hxap  hxcabg  hxcan  hxdiam  hxhyp  hxicd  hxpcor  hxstr hxpmi  hxpaim smoke dbp sbp hr  bminum ef age sex  nyha; 
		run;
		proc phreg data = itt; 
			model t2dthc*death(0) = rand /rl; 
			title 'ITT, no adjustment (Cox model)';
		run;

		%if &row = 1 %then %do;
			proc phreg data = itt; 
				model t2dthc*death(0) = rand /rl; 
				strata study_;
				title 'ITT, substudy only (row 1, Cox model)';
				ods output parameterestimates = param1;
			run;

			proc printto print = "Table1_ITT.Cox.&nboot.boot.rtf";
			run;

			data param1; 
				set param1; 		
				model_type = "cox"; 
			run;

			proc print data = param1;
				var Parameter HazardRatio HRLowerCL HRUpperCL ProbChiSq model_type; 
				title 'ITT, substudy only (row 1, Cox model)';
			run;

			proc printto;
			run;
		%end;
		%else %if &row = 2 %then %do;
			proc phreg data = itt; 
				class nyha (ref = "4") / param = ref;
				model t2dthc*death(0) = rand dgbb_b dgdiu_b dgacei_b dglild_b	
					hxaf  hxap  hxcabg  hxcan  hxdiam  hxhyp  hxicd  hxpcor hxstr  hxpmi  hxpaim smoke dbp sbp hr  bminum ef age sex  nyha /rl; 
				strata study_;
				title 'ITT, substudy and covariates (row 2, Cox model)';
				ods output parameterestimates = param2;
			run;

			proc printto print = "Table1_ITT.Cox.&nboot.boot.rtf";
			run;

			data pvalues1;
			set param2;
				if Parameter in ("age0", "hr0", "dbp0", "sbp0") then do; 
					 HazardRatio_new = exp(Estimate*10);
					 HRLowerCL_new = exp(Estimate*10 - 1.96*StdErr);
			 		 HRUpperCL_new = exp(Estimate*10 + 1.96*StdErr);
				end;
				else do;
					HazardRatio_new = HazardRatio;
					HRLowerCL_new = HRLowerCL;
				 	HRUpperCL_new = HRUpperCL;
				end;
				model_type = "cox";
				keep Parameter HazardRatio_new HRLowerCL_new HRUpperCL_new ProbChiSq StdErr model_type ;
			run;

			proc print data = pvalues1;
				var Parameter HazardRatio_new HRLowerCL_new HRUpperCL_new model_type; 
				title 'ITT, substudy and covariates (row 2, Cox model)';
			run;

			proc printto;
			run;
		%end;
	%end;
	%else %if &row ge 3 %then %do;
		/*Pooled logistic regression, hazard ratios only*/

		%let outdest = "Table1_ITT.PooledLogistic.&nboot.boot.rtf";
		%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compl_t_c3, varout_baseline = compl_t0,
			id = b_patid , time = , completecase = 0, compyes =0,dichot = 2, outdatatype = 3, 
			LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

		%if &row = 3 %then %do;
			%stdz(outdest = &outdest, inset =cuminc,  title_main = "ITT, no adjustment (row 3, Pooled logistic)", 
				nboot= &nboot,  expvar0 = compl_t0, expvar = rand,  arm = ., eof = &eof, pha = 1,
				event = event, measurevar = comp_m, censLTFU =  censLTFU3, ipcw = 0, iptw = 0, 
				timevar= weektime, visitvar = visit_t , recode = ., adjust = 0, itt = 1, ppe = 0, crude = 1, protocol = 0, cuminc = 0);
		%end;
		%else %if &row = 4 %then %do;
			%stdz(outdest = &outdest, inset =cuminc,  title_main = "ITT, substudy only (row 4, Pooled logistic)", 
				nboot= &nboot,  expvar0 = compl_t0, expvar = rand,  arm = ., eof = &eof,pha = 1,
				event = event, measurevar = comp_m, censLTFU = ., ipcw = 0, iptw = 0, 
				timevar= weektime, visitvar = visit_t , recode = ., adjust = 0, itt = 1, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
				covs = study_);
		%end;
		%else %if &row = 5 %then %do;
			%stdz(outdest = &outdest, inset =cuminc,  title_main = "ITT, substudy and covariates (row 5, Pooled logistic)", 
				nboot= &nboot,  expvar0 = compl_t0, expvar = rand,  arm = ., eof = &eof,pha = 1,
				event = event, measurevar = comp_m, censLTFU =  censLTFU3, ipcw = 0, iptw = 0, 
				timevar= weektime, visitvar = visit_t , recode = ., adjust = 0, itt = 1, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
				covs = &covs2);
		%end;
	%end;
%mend Table1;	


/**************************************************/
/*Macro Table 2: Placebo Arm***********************/
/**************************************************/

%macro Table2(nboot = , eof = , row = , sens = 0, sens_option = );
  	%if &row = 1 %then %do;
		/*Row 1: Replication of Granger et al, as described in methods*/

		%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_m0, varout_baseline = compyes_t0,
			id = b_patid , time =weektime , completecase = 0, compyes = 1, dichot = 0.8, outdatatype = 2, 
			recode = 1, recode_var1 = comp_m, recodeifeq_cond1 = 0, recodeto_cond1 = 3, norecodebefore = 3);

		%if &sens = 0 %then %do;
			%let outdest = "Table2_row1.MainAnalysis.&nboot.boot.rtf";
			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 1: Pooled logistic: Replication of methods, with binary adherence", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes, arm = 0,  eof = &eof, pha = 1,
				event = event_new, ipcw = 0, iptw =0, timevar= weektime, visitvar = visit_t , recode = -3, 
				adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, covs =  &covs1, nocens = 1, outc_only = 1);
		%end;
		%else %if &sens = 1 %then %do;
			%let outdest = "Table2_row1.SensitivityAnalyses.&nboot.boot.rtf";		
			%if &sens_option = "Model type" %then %do;
				%nostdz(inset = cuminc, outdest = &outdest, title_main = 'Cox: Replication of methods, with binary adherence',  expvar =compl, crude = 0 ,
				  	class= 0, id = b_patid, logit = 0,  event =event_new, itt = 0, ppe = 0, subset = rand, restrict_subset = 0, covs = &covs1, 
					recode = 1 , stop = 1, start = 0, pooled_l = 0, ag = 1);
			%end;

			%if &sens_option = "Adherence missingness" %then %do;
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_m0, varout_baseline = compyes_t0,
					id = b_patid , time =weektime , completecase = 0, compyes = 1, dichot = 0.8, outdatatype = 2, 
					recode = 1, recode_var1 = comp_m, recodeifeq_cond1 = 0, recodeto_cond1 = 3, norecodebefore = 3);

				%nostdz(inset = cuminc, outdest = &outdest, title_main = 'Cox (Closest replication): drop if miss, model linear',  expvar =compl, crude = 0 ,
					class= 0, id = b_patid, logit = 0,  event =event_new, itt = 0, ppe = 0, subset = rand, restrict_subset = 0, covs =  &covs1, 
					recode = 2 , stop = 0, start = 0, pooled_l = 0, ag = 1);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Pooled logistic: drop if miss, model linear", 
					nboot= &nboot,  expvar0 = compl_t0, expvar = compl, arm = 0,  eof = &eof, pha = 1,
					event = event_new, ipcw = 0, iptw =0, timevar= weektime, visitvar = visit_t , recode = -2, 
					adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, covs =  &covs1, nocens = 1, outc_only = 1);
			%end;
		%end;
	%end;
   	%else %if &row = 2 %then %do;
		/*Row 2: Improved covariates*/

		%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_m0, varout_baseline = compyes_t0,
			id = b_patid , time =weektime , completecase = 0, compyes = 1, dichot = 0.8, outdatatype = 2, 
			recode = 1, recode_var1 = comp_m, recodeifeq_cond1 = 0, recodeto_cond1 = 3, norecodebefore = 3);

		%if &sens = 0 %then %do;
			%let outdest = "Table2_row2.MainAnalysis.&nboot.boot.rtf";

			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 2: Pooled logistic: Improved covariate adjustment", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes, arm = 0,  eof = &eof, pha = 1,
				event = event_new, ipcw = 0, iptw =0, timevar= weektime, visitvar = visit_t , recode = -3, 
				adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, covs =  &covs2, nocens = 1, outc_only = 1);
		%end;
		%else %if &sens = 1 %then %do;
			%let outdest = "Table2_row2.SensitivityAnalyses.&nboot.boot.rtf";

			%if &sens_option = "Model type" %then %do;
				%nostdz(inset = cuminc, outdest = &outdest, title_main = 'Cox: Replication of methods, with binary adherence',  expvar =compl, crude = 0 ,
				  	class= 2, id = b_patid, logit = 0,  event =event_new, itt = 0, ppe = 0, subset = rand, restrict_subset = 0, covs = &covs2, 
					recode = 1 , stop = 1, start = 0, pooled_l = 0, ag = 1);
			%end;
		%end;
	%end;
	%else %if &row = 3 %then %do;
		/*Row 3: Pooled logistic regression with improved covariate modeling & standardization*/

		%if &sens = 0 %then %do;
			%let outdest = "Table2_row3.MainAnalysis.&nboot.boot.rtf";
	
			%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_m0, varout_baseline = compyes_t0,
				id = b_patid , time =weektime , completecase = 0, compyes = 1, dichot = 0.8, outdatatype =2, 
				recode = 1, recode_var1 = comp_m, recodeifeq_cond1 = 0, recodeto_cond1 = 3, norecodebefore = 3);
	
			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 3: with standardization", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes,  arm = 0,  eof = &eof, pha = 1,
				event = event_new, ipcw = 0, iptw =0, timevar= weektime, visitvar = visit_t , recode = -3, 
				adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, covs =  &covs2, nocens = 1);
		%end;
	%end;

	%else %if &row = 4 %then %do;
		/*Row 4: Censor when lost to follow-up -- 3 consecutive missed visits*/
		%if &sens = 0 %then %do;
			%let outdest = "Table2_row4.MainAnalysis.&nboot.boot.rtf";

			%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
				id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
				LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: censor when lost to follow-up (>3 missed visits)", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0,  eof = &eof, pha = 1,
				event = event, measurevar = comp_m, censLTFU = censLTFU3, ipcw = 0, iptw =0, 
				timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
				covs = &covs2, nocens = 1);
		%end;
		%else %if &sens = 1 %then %do;
			%let outdest = "Table2_row4.SensitivityAnalyses.&nboot.boot.rtf";
			%if &sens_option = "Censoring time" %then %do;
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: censor when lost to follow-up (>1 missed visits)", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0, eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU1, ipcw = 0, iptw =0, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs = &covs2, nocens = 1);

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: censor when lost to follow-up (>2 missed visits)", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0, eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU2, ipcw = 0, iptw =0, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs = &covs2, nocens = 1);

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: censor when lost to follow-up (>4 missed visits)", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0, eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU4, ipcw = 0, iptw =0, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs = &covs2, nocens = 1);
			%end;
		%end;
	%end;
 	%if &row = 5 %then %do; /*Additionally censor when adherence changes from baseline*/
		%if &sens = 0 %then %do;
			%let outdest = "Table2_row5.MainAnalysis.&nboot.boot.rtf";

			%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
				id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
				LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 5: Censor when adherence changes", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0,  eof = &eof, pha = 1,
				event = event, measurevar = comp_m, censLTFU = censLTFU3, ipcw = 0, iptw =0, 
				timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
				covs =  &covs2);
		%end;
		%else %if &sens = 1 %then %do;
			%let outdest = "Table2_row5.SensitivityAnalyses.&nboot.boot.rtf";

			%if &sens_option = "Censoring time" %then %do;

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 5: Censor > 1 missed visit", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU1, ipcw = 0, iptw =0, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);
	
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 5: Censor > 2 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU1, ipcw = 0, iptw =0, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);
	
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 5: Censor > 4 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU1, ipcw = 0, iptw =0, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);
			%end;
		%end;
   	%end;
   	%if &row = 6 %then %do;
		%if &sens = 0 %then %do;
			%let outdest = "Table2_row6.MainAnalysis.&nboot.boot.rtf";

			%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
				id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
				LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 6: Adjustment for postrandomization covariates", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3, arm = 0,  eof = &eof, pha = 1,
				event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 1, iptw = 1, 
				timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
				covs =  &covs2, covs_t = &covs_t);
		%end;
		%else %if &sens = 1 %then %do;
			%let outdest = "Table2_row6.SensitivityAnalyses.&nboot.boot.rtf";

			%if &sens_option = "Censoring time" %then %do;

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 6: Censor > 1 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3, arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU1, n_cens = 1, ipcw = 1, iptw = 1, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);
		
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 6: Censor > 2 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3, arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU2, n_cens = 2, ipcw = 1, iptw = 1, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 6: Censor > 4 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3, arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU4, n_cens = 4, ipcw = 1, iptw = 1, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);

			%end;
			%else %if &sens_option = "IP weights" %then %do;

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 6: IPTW only", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3, arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 0, iptw = 1, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);
	
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 6: IPCW only", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3, arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 1, iptw = 0, 
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);
			%end;
		%end;
  	%end;
   	%if &row = 7 %then %do;
		%if &sens = 0 %then %do;
			%let outdest = "Table2_row7.MainAnalysis.&nboot.boot.rtf";

			%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
				id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
				LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: Dose-response, with baseline standardization", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1, 
				event = event, measurevar = comp_m, censLTFU = censLTFU3, ipcw = 0, iptw =0, contA = 1, doseform = 4,
				timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
				covs =  &covs2);
		%end;
		%else %if &sens = 1 %then %do;
			%let outdest = "Table2_row7.SensitivityAnalyses.&nboot.boot.rtf";

			%if &sens_option = "Censoring time" %then %do;

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: Censor > 1 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU1, ipcw = 0, iptw =0, contA = 1, doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);
	
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: Censor > 2 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU2, ipcw = 0, iptw =0, contA = 1, doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: Censor > 4 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU4, ipcw = 0, iptw =0, contA = 1, doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);
			%end;
			%else %if &sens_option = "Functional form" %then %do;

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
			 		LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: At+c(At-1)", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, ipcw = 0, iptw =0, contA = 1, doseform = 3,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: c(At)", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, ipcw = 0, iptw =0, contA = 1, doseform = 1,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: c(At) + c(At)^2 ", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, ipcw = 0, iptw =0, contA = 1, doseform = 2,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2);
			%end;
		%end;
	%end;
   	%if &row = 8 %then %do;
		%if &sens = 0 %then %do;
			%let outdest = "Table2_row8.MainAnalysis.&nboot.boot.rtf";
			%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
				id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
				LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

			%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row8: Dose-response, with baseline standardization", 
				nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
				event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 1, iptw = 1,  contA = 1, doseform = 4,
				timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
				covs =  &covs2, covs_t = &covs_t);
		%end;
		%else %if &sens = 1 %then %do;
			%let outdest = "Table2_row8.SensitivityAnalyses.&nboot.boot.rtf";
			
			%if &sens_option = "Censoring time" %then %do;
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 8: Censor > 1 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU1, n_cens = 1, ipcw = 1, iptw = 1,  contA = 1,  doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);
	
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 8: Censor > 2 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU2, n_cens = 2, ipcw = 1, iptw = 1,  contA = 1,  doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);
	
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 8: Censor > 4 missed visits", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU4, n_cens = 4, ipcw = 1, iptw = 1,  contA = 1,  doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);
			%end;
			%else %if &sens_option = "IP weights" %then %do;

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
			 		LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row8: IPTW only", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 0, iptw = 1,  contA = 1,  doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row8: IPCW only", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 1, iptw = 0,  contA = 1,  doseform = 4,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);	
			%end;
			%else %if &sens_option = "Functional form" %then %do;

				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
		 			LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 8: At+c(At-1)", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 1, iptw = 1,  contA = 1, doseform = 3,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 8: c(At)", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 1, iptw = 1,  contA = 1, doseform = 1,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);

				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 8: c(At) + c(At)^2", 
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3_avg,  arm = 0,  eof = &eof, pha = 1,
					event = event, measurevar = comp_m, censLTFU = censLTFU3, n_cens = 3, ipcw = 1, iptw = 1,  contA = 1, doseform = 2,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
					covs =  &covs2, covs_t = &covs_t);
			%end;
		%end;
   	%end;
%mend Table2;

/**************************************************/
/*Macro Table 3: Per-protocol Effects**************/
/**************************************************/

%macro Table3(nboot = , eof = , effect = , adjust = , est = , sens = 0, sens_option = , graph = 0);
	%if &est = "HR" %then %do;
		%let pha = 1;
	%end;
	%else %if &est = "RD" %then %do;
		%let pha = 0;
	%end;

	%if &effect = "ITT" %then %do;
		%if &graph = 0 %then %do;
			%let outdest = "Table3_row0.MainAnalysis.&nboot.boot.rtf";
		%end;
		%else %if &graph = 1 %then %do;
			%let outdest = "Figure1.ITT.&nboot.boot.rtf";
		%end;
		%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compl_t_c3, varout_baseline = compl_t0,
			id = b_patid , completecase = 0, compyes =0,dichot = 2, outdatatype = 3, 
			LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

		%stdz(outdest = &outdest, inset =cuminc,  title_main = "ITT, substudy and covariates (row 0)", 
			nboot= &nboot,  expvar0 = compl_t0, expvar = rand, eof = &eof, pha = &pha,
			event = event, measurevar = comp_m, censLTFU =  censLTFU3, ipcw = 0, iptw = 0, graph = &graph,
			timevar= weektime, visitvar = visit_t , adjust = 0, itt = 1, ppe = 0, crude = 0, protocol = 0, cuminc = 0, 
			covs = &covs2);
	%end;
		
	%if &effect = "PPE, major AE" %then %do;
		%if &sens = 0 %then %do;
			%if &adjust = "dose-response" %then %do;
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = weektime, completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
					recodeifeq_cond2 = 2, recodeto_cond2 = 0,  recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
					LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);
			%end;
			%else %do;
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
					n_cens_cond = 2, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 2, 
					censif_var2 = everdisc_r, censif_cond2 = 3, nocensbefore = 0,
					LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);
			%end;
			%if &adjust = "none" %then %do;
				%let outdest = "Table3_row1.MainAnalysis.&nboot.boot.rtf";
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 1: PPE, major AEs, unadjusted", 
					nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha = &pha,
					event = event, measurevar = comp_m, censLTFU =  censLTFU3, ipcw = 0, iptw = 0, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 1, crude = 1, protocol = 1, cuminc = 0);
			%end;
			%else %if &adjust = "substudy" %then %do;
				%let outdest = "Table3_row2.MainAnalysis.&nboot.boot.rtf";
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 2: PPE, major AEs, substudy adjusted", eof = &eof, pha = &pha,
					nboot= &nboot,  expvar0 = compyes_t0, expvar =compyes_t_c3,  arm = ., 
					event = event, measurevar = comp_m, censLTFU = ., ipcw = 0, iptw = 0, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 1, crude = 0, protocol =1, cuminc = 0, 
					covs = study_);
			%end;
			%else %if &adjust = "covs0" %then %do;
				%let outdest = "Table3_row3.MainAnalysis.&nboot.boot.rtf";
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 3: PPE, major AEs, baseline adjusted", eof = &eof, pha = &pha,
					nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3, arm = ., 
					event = event, measurevar = comp_m, censLTFU =  censLTFU3, ipcw = 0, iptw = 0, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 1, crude = 0, protocol = 1, cuminc = 0,
					covs = &covs2);

			%end;
			%else %if &adjust = "covs_t" %then %do;
				%if &graph = 0 %then %do;
					%let outdest = "Table3_row4.MainAnalysis.&nboot.boot.rtf";
				%end;
				%else %if &graph = 1 %then %do;
					%let outdest = "Figure1.PPEmajor.&nboot.boot.rtf";
				%end;
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: PPE, major AEs, post-baseline adjusted", eof = &eof, pha = &pha,
					nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., graph = &graph,
					event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 1, cuminc = 0, 
					covs = &covs2, covs_t = &covs_t);
			%end;
			%else %if &adjust = "dose-response" %then %do;
				%let outdest = "Table3_row5.MainAnalysis.&nboot.boot.rtf";
				%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: PPE, major AEs, dose-response", 
					nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 0, iptw = 1, protocol = 1, doseform = 4,
					timevar= weektime, visitvar = visit_t ,  eof =&eof, pha = &pha,
					covs = &covs2, covs_t = &covs_t);
			%end;
		%end;
		%if &sens = 1 %then %do;
			%if &sens_option = "Censoring time" %then %do;				
				%if &adjust = "covs_t" %then %do;
					%let outdest = "Table3_row4.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 2, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 2, 
						censif_var2 = everdisc_r, censif_cond2 = 3, nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);
					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: Censor > 1 missed visit", eof = &eof, pha = &pha,
						nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., 
						event = event, measurevar = comp_m, censLTFU = censLTFU1,  n_cens = 1, ipcw = 1, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 1, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 2, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 2, 
						censif_var2 = everdisc_r, censif_cond2 = 3, nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);
					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: Censor > 2 missed visits", eof = &eof, pha = &pha,
						nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., 
						event = event, measurevar = comp_m, censLTFU = censLTFU2,  n_cens = 2, ipcw = 1, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 1, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1, dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 2, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 2, 
						censif_var2 = everdisc_r, censif_cond2 = 3, nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);
					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: Censor > 4 missed visits", eof = &eof, pha = &pha,
						nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., 
						event = event, measurevar = comp_m, censLTFU = censLTFU4,  n_cens = 4, ipcw = 1, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 1, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);
				%end;
				%else %if &adjust = "dose-response" %then %do;
					%let outdest = "Table3_row5.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime, completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 0,  recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);
					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: Censor > 1 missed visit", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU1,  n_cens = 1, ipcw = 1, iptw = 1, protocol = 1, doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime, completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 0,  recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);
					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: Censor > 2 missed visits", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU2,  n_cens = 2, ipcw = 1, iptw = 1, protocol = 1, doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime, completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 0,  recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);
					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: Censor > 4 missed visits", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU4,  n_cens = 4, ipcw = 1, iptw = 1, protocol = 1, doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);
				%end;
			%end;
			%else %if &sens_option = "IP weights" %then %do;
				%if &adjust = "covs_t" %then %do;
					%let outdest = "Table3_row4.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 2, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 2, 
						censif_var2 = everdisc_r, censif_cond2 = 3, nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: IPTW only", eof = &eof, pha = &pha,
						nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., 
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 0, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 1, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);

					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 4: IPCW only", eof = &eof, pha = &pha,
						nboot= &nboot,  expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., 
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 0, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 1, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);
				%end;
				%else %if &adjust = "dose-response" %then %do;
					%let outdest = "Table3_row5.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime, completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 0,  recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);
	
					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: IPTW only", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 0, iptw = 1, protocol = 1,doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: IPCW only", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 0, protocol = 1,doseform = 4,
						timevar= weektime, visitvar = visit_t , eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);
				%end;
			%end;
			%else %if &sens_option = "Functional form" %then %do;
				%if &adjust = "dose-response" %then %do;
					%let outdest = "Table3_row5.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime, completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 0,  recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: c(At)", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, protocol = 1, doseform = 1,
						timevar= weektime, visitvar = visit_t ,  eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: c(At) + c(At)2", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, protocol = 1, doseform = 2,
						timevar= weektime, visitvar = visit_t ,  eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 5: At + c(At-1)", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, protocol = 1, doseform = 3,
						timevar= weektime, visitvar = visit_t , eof =&eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);
				%end;
			%end;
		%end;
	%end;
	%else %if &effect = "PPE, any AE" %then %do;
		%if &sens = 0 %then %do;
			%if &adjust = "dose-response" %then %do;
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = weektime , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
					recodeifeq_cond2 = 2, recodeto_cond2 = 1, recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
					LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);
			%end;
			%else %do;
				%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
					id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
					n_cens_cond = 1, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 3,  nocensbefore = 0,
					LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);
			%end;
			
			%if &adjust = "none" %then %do;
				%let outdest = "Table3_row6.MainAnalysis.&nboot.boot.rtf";
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 6: PPE, any AEs, unadjusted", 
					nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha = &pha,
					event = event, measurevar = comp_m, censLTFU =  censLTFU3, ipcw = 0, iptw = 0, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 1, crude = 1, protocol = 2, cuminc = 0);
			%end;
			%else %if &adjust = "substudy" %then %do;
				%let outdest = "Table3_row7.MainAnalysis.&nboot.boot.rtf";
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 7: PPE, any AEs, substudy adjusted", 
					nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha = &pha,
					event = event, measurevar = comp_m, censLTFU = ., ipcw = 0, iptw = 0, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 1, crude = 0, protocol =2, cuminc = 0, 
				covs = study_);
			%end;
			%else %if &adjust = "covs0" %then %do;
				%let outdest = "Table3_row8.MainAnalysis.&nboot.boot.rtf";
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 8: PPE, any AEs, baseline adjusted", 
					nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha = &pha, 
					event = event, measurevar = comp_m, censLTFU =  censLTFU3, ipcw = 0, iptw = 0, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 0, itt = 0, ppe = 1, crude = 0, protocol = 2, cuminc = 0,
					covs = &covs2);
			%end;
			%else %if &adjust = "covs_t" %then %do;
				%if &graph = 0 %then %do;
					%let outdest = "Table3_row9.MainAnalysis.&nboot.boot.rtf";
				%end;
				%else %if &graph = 1 %then %do;
					%let outdest = "Figure1.PPEany.&nboot.boot.rtf";
				%end;
				%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 9: PPE, any AEs, post-baseline adjusted", 
					nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha = &pha, graph = &graph,
					event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, censvarname = cens_prot,
					timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 2, cuminc = 0, 
					covs = &covs2, covs_t = &covs_t);
			%end;
			%else %if &adjust = "dose-response" %then %do;
				%let outdest = "Table3_row10.MainAnalysis.&nboot.boot.rtf";
				%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: PPE, any AEs, dose-response", 
					nboot= &nboot, expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
				  	event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, protocol = 2, doseform = 4,
					timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
					covs = &covs2, covs_t = &covs_t);
			%end;
		%end;		
		%else %if &sens = 1 %then %do;
			%if &sens_option = "Censoring time" %then %do;				
				%if &adjust = "covs_t" %then %do;
					%let outdest = "Table3_row9.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 1, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 3,  nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);
					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 9: Censor > 1 missed visit", 
						nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha =  &pha,
						event = event, measurevar = comp_m, censLTFU = censLTFU1,  n_cens = 1, ipcw = 1, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 2, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 1, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 3,  nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);
					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 9: Censor > 2 missed visits", 
						nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha =  &pha,
						event = event, measurevar = comp_m, censLTFU = censLTFU2,  n_cens = 2, ipcw = 1, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 2, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 1, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 3,  nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);
					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 9: Censor > 4 missed visits", 
						nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha =  &pha,
						event = event, measurevar = comp_m, censLTFU = censLTFU4,  n_cens = 4, ipcw = 1, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 2, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);
				%end;
				%else %if &adjust = "dose-response" %then %do;
					%let outdest = "Table3_row10.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 1, recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU1, LTFUaftervisit = 1, LTFUcond_var1 = vs_adherence);
					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: Censor > 1 missed visit", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU1,  n_cens = 1, ipcw = 1, iptw = 1, protocol = 2, doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 1, recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU2, LTFUaftervisit = 2, LTFUcond_var1 = vs_adherence);
					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: Censor > 2 missed visits", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU2,  n_cens = 2, ipcw = 1, iptw = 1, protocol = 2, doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 1, recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU4, LTFUaftervisit = 4, LTFUcond_var1 = vs_adherence);
					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: Censor > 4 missed visits", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU4,  n_cens = 4, ipcw = 1, iptw = 1, protocol = 2, doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);
				%end;
			%end;
			%else %if &sens_option = "IP weights" %then %do;
				%if &adjust = "covs_t" %then %do;
					%let outdest = "Table3_row9.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						n_cens_cond = 1, censvarname = cens_prot, censif_var1 = everdisc_r, censif_cond1 = 3,  nocensbefore = 0,
						LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 9: IPTW only", 
						nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha = &pha,
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw =0, iptw = 1, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 2, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);

					%stdz(outdest = &outdest, inset =cuminc,  title_main = "Row 9: IPCW only", 
						nboot= &nboot, expvar0 = compyes_t0, expvar = compyes_t_c3,  arm = ., eof = &eof, pha = &pha,
						event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw =1, iptw = 0, censvarname = cens_prot,
						timevar= weektime, visitvar = visit_t , recode = 2, adjust = 1, itt = 0, ppe = 1, crude = 0, protocol = 2, cuminc = 0, 
						covs = &covs2, covs_t = &covs_t);
				%end;
				%else %if &adjust = "dose-response" %then %do;
					%let outdest = "Table3_row10.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 1, recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: IPTW only", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 0, iptw = 1, protocol = 2,doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: IPCW only", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 0, protocol = 2,doseform = 4,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);
				%end;
			%end;
			%else %if &sens_option = "Functional form" %then %do;
				%if &adjust = "dose-response" %then %do;
					%let outdest = "Table3_row10.SensitivityAnalyses.&nboot.boot.rtf";
					%cumavg(inset = charm.final_weeks1, outset = cuminc, varname_out = compyes_t_c3, varout_baseline = compyes_t0,
						id = b_patid , time = weektime , completecase = 0, compyes =1,dichot = 0.8, outdatatype = 3, 
						recode = 2, recode_var1 = everdisc_r, recodeifeq_cond1 = 1, recodeto_cond1 = 1, 
						recodeifeq_cond2 = 2, recodeto_cond2 = 1, recodeifeq_cond3 = 3, recodeto_cond3 = 0, norecodebefore = 0,
						LTFU = 1, LTFUvarname = censLTFU3, LTFUaftervisit = 3, LTFUcond_var1 = vs_adherence);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: c(At)", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, protocol = 2,doseform = 1,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: c(At) + c(At)2", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, protocol = 2,doseform = 2,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);

					%dose_r(outdest = &outdest, inset =cuminc,  title_main = "Row 10: At + c(At-1)", 
						nboot= &nboot,   expvar0 = compyes_t0, expvar_eof = compyes_t_c3_avg_eof, expvar = compyes_t_c3_avg,    
					  	event = event, measurevar = comp_m, censLTFU = censLTFU3,  n_cens = 3, ipcw = 1, iptw = 1, protocol = 2,doseform = 3,
						timevar= weektime, visitvar = visit_t ,  eof = &eof, pha = &pha,
						covs = &covs2, covs_t = &covs_t);
				%end;
			%end;
		%end;
	%end;
%mend Table3;




