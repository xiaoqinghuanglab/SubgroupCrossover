options validvarname=v7 mprint mlogic symbolgen nodate nonumber;
options compress=yes reuse=yes msglevel=i;

libname NACC "/home/u63719145/Clinical Paper/SubgroupEffect/real data";

proc datasets lib=work nolist kill; quit;

/*-------------------------------*
 | Build COHORT from SW_NOAD_CLEAN
 *-------------------------------*/
data _with_state0;
    set nacc.M_COMBM_CLEAN;

    length dname $300 state4 $12;
    f_DONEPEZIL = 0;
    f_DONE_MEM  = 0;
    f_MEMANTINE = 0;

    array drugs[*] $ DRUG1-DRUG40;
    do _k = 1 to dim(drugs);
        dname = upcase(compbl(strip(drugs[_k])));
        if dname ne "" then do;
            if prxmatch('/DONEPEZIL.*MEMANTINE|MEMANTINE.*DONEPEZIL|NAMZARIC/', dname) then f_DONE_MEM = 1;
            else do;
                if prxmatch('/\b(DONEPEZIL|ARICEPT)\b/', dname) then f_DONEPEZIL = 1;
                if prxmatch('/\b(MEMANTINE|NAMENDA)\b/', dname) then f_MEMANTINE = 1;
            end;
        end;
    end;

    if (f_DONE_MEM = 1 or (f_DONEPEZIL = 1 and f_MEMANTINE = 1)) then state4 = "DONEP_MEM";
    else if (f_MEMANTINE = 1 and f_DONEPEZIL = 0 and f_DONE_MEM = 0) then state4 = "MEM_ONLY";
    else state4 = "NO_AD";

    drop _k dname f_DONEPEZIL f_DONE_MEM f_MEMANTINE;
run;

/* one row per subject-phase:
   prefer last non-missing CDRSUM row within each phase */
data _with_state0;
    set _with_state0;
    _hasy = (not missing(CDRSUM));
run;

proc sort data=_with_state0 out=_with_state;
    by NACCID phase descending _hasy visit_date;
run;

data _with_state;
    set _with_state;
    by NACCID phase;
    if first.phase;
    drop _hasy;
run;

/* keep complete 3-phase subjects only */
proc sql;
    create table _okid as
    select NACCID
    from _with_state
    group by NACCID
    having sum(phase=1) > 0 and sum(phase=2) > 0 and sum(phase=3) > 0;
quit;

proc sql noprint;
    create table _dir as
    select NACCID,
           case when state4 = "MEM_ONLY" then "M_to_C" else "C_to_M" end as direction length=6,
           case when state4 = "MEM_ONLY" then "M_to_C" else "C_to_M" end as switch_seq length=6
    from _with_state
    where phase = 2;
quit;

proc sort data=_with_state; by NACCID; run;
proc sort data=_dir; by NACCID; run;
proc sort data=_okid; by NACCID; run;

data cohort;
    merge _with_state(in=a) _dir(in=b) _okid(in=c);
    by NACCID;
    if a and b and c;

    tripleid  = NACCID;
    visitnum3 = phase - 1;
    BASELINE = (phase = 1);
    PERIOD2  = (phase = 3);
    analysis_y = CDRSUM;

    length current_trt $8;
    select (phase);
        when (1) do; ACTIVE = 0; CARRYOVER = 0; current_trt = "NONE"; end;
        when (2) do;
            CARRYOVER = 0;
            if state4 = "MEM_ONLY" then do; ACTIVE = 0; current_trt = "MEM"; end;
            else do; ACTIVE = 1; current_trt = "COMBO"; end;
        end;
        when (3) do;
            if state4 = "MEM_ONLY" then do; ACTIVE = 0; CARRYOVER = 1; current_trt = "MEM"; end;
            else do; ACTIVE = 1; CARRYOVER = 0; current_trt = "COMBO"; end;
        end;
        otherwise;
    end;
run;


title "Cohort QC (expect ~3 rows per subject)";
proc sql;
    select count(distinct NACCID) as n_subjects format=comma12.,
           count(*) as n_rows format=comma12.
    from cohort;
quit;
title;

proc datasets lib=work nolist; delete _with_state0 _with_state _dir _okid; quit;

/*-------------------------------*
 | Helper macros
 *-------------------------------*/
%macro tag_append(ds=, out=, var=, model=, mean=, src=);
    %if %sysfunc(exist(&ds)) %then %do;
        data _tag;
            length subg_var $32 model $12 mean_model $24 data_src $12
                   Effect $32 Label $80 Descr $80;
            set &ds;
            subg_var="%upcase(&var)"; model="&model"; mean_model="&mean"; data_src="&src";
        run;
        proc append base=&out data=_tag force; run;
        proc datasets lib=work nolist; delete _tag; quit;
    %end;
%mend;

%macro get_ll(name=);
    %if %sysfunc(exist(mlfit_&name)) %then %do;
        data ll_&name;
            set mlfit_&name;
            if index(Descr,"-2 Log Likelihood")>0 then do; neg2ll=Value; output; end;
            keep neg2ll;
        run;
    %end;
    %else %do;
        data ll_&name; neg2ll=.; output; run;
    %end;
%mend;

%macro make_lrt(var=, src=, cov=, test=, full=, reduced=, df=);
    data lrt_one;
        length subg_var $32 data_src $12 Covariance $12 Test $80;
        merge ll_&full(rename=(neg2ll=neg2ll_full))
              ll_&reduced(rename=(neg2ll=neg2ll_reduced));
        subg_var="%upcase(&var)"; data_src="&src"; Covariance="&cov"; Test="&test"; DF=&df;
        LR_Stat=neg2ll_reduced-neg2ll_full;
        invalid_lrt=0;
        if not missing(LR_Stat) and LR_Stat<0 then do;
            if abs(LR_Stat)<1e-6 then LR_Stat=0; else invalid_lrt=1;
        end;
        if invalid_lrt=0 and not missing(LR_Stat) then P_Value=1-probchi(LR_Stat,DF);
        else P_Value=.;
    run;
    proc append base=lrt_all data=lrt_one force; run;
%mend;

%macro fit_mcombm(var=, pos=1, neg=0, sgrow=BASE);
    %local src_label unit tvar idx dsid rc reqvars missing_req _i _v n_unit n_subj n_seq minc n_cell;
    %let src_label=M_COMBM;
    %let var=%upcase(&var);
    %let sgrow=%upcase(&sgrow);
    %let unit=tripleid;
    %let tvar=visitnum3;
    %if &sgrow=BASE %then %let idx=(visitnum3=0);
    %else %if &sgrow=P1 %then %let idx=(visitnum3=1);
    %else %do; %put ERROR: sgrow must be BASE/P1.; %return; %end;

    %let reqvars=NACCID tripleid visitnum3 switch_seq direction current_trt BASELINE PERIOD2 ACTIVE CARRYOVER analysis_y &var;
    %let missing_req=0;
    %let dsid=%sysfunc(open(cohort,i));
    %if &dsid=0 %then %do; %put ERROR: cannot open COHORT.; %return; %end;
    %do _i=1 %to %sysfunc(countw(&reqvars));
        %let _v=%scan(&reqvars,&_i);
        %if %sysfunc(varnum(&dsid,&_v))=0 %then %do;
            %put ERROR: missing variable &_v in COHORT.;
            %let missing_req=1;
        %end;
    %end;
    %let rc=%sysfunc(close(&dsid));
    %if &missing_req=1 %then %return;

    proc datasets lib=work nolist;
        delete ana_long unit_level subg_src subg_unit raw_freq screen_one cell_cnt _cell_raw _seq_levels _subg_levels _cell_grid
               fit_all est_all typ_all cov_all lrt_all key_est key_typ
               fit_cs est_cs typ_cs cov_cs fit_un est_un typ_un cov_un fit_csn est_csn typ_csn cov_csn fit_unn est_unn typ_unn cov_unn
               mlfit_: ll_: lrt_one;
    quit;

    data subg_src; set cohort; if &idx; subg_pre=&var; keep &unit subg_pre; run;
    proc sort data=subg_src nodupkey; by &unit; run;

    proc freq data=subg_src noprint;
        tables subg_pre / missing out=raw_freq;
    run;

    data raw_freq;
        length subg_var $32 data_src $12 sgrow $8;
        set raw_freq;
        subg_var = "&var";
        data_src = "&src_label";
        sgrow = "&sgrow";
    run;

    data subg_unit;
        set subg_src;
        if subg_pre in (&neg) then SUBG=0;
        else if subg_pre in (&pos) then SUBG=1;
        else delete;
    run;
    proc sort data=subg_unit; by &unit; run;
    proc sort data=cohort; by &unit; run;

    data ana_long;
        merge cohort(in=a) subg_unit(in=b);
        by &unit;
        if a and b;
        ACTIVE_SUBG=ACTIVE*SUBG;
        CARRYOVER_SUBG=CARRYOVER*SUBG;
    run;

    proc sort data=ana_long(where=&idx) out=unit_level nodupkey; by &unit; run;
    proc sql noprint;
        select count(distinct &unit), count(distinct NACCID), count(distinct switch_seq)
            into :n_unit, :n_subj, :n_seq from unit_level;
        create table _cell_raw as
            select switch_seq, SUBG, count(distinct &unit) as n_units
            from unit_level group by switch_seq,SUBG;

        create table _seq_levels as
            select distinct switch_seq
            from unit_level
            where not missing(switch_seq);
    quit;

    data _subg_levels;
        SUBG = 0; output;
        SUBG = 1; output;
    run;

    proc sql;
        create table _cell_grid as
        select s.switch_seq, g.SUBG
        from _seq_levels as s, _subg_levels as g;

        create table cell_cnt as
        select
            "&src_label" as data_src length=12,
            "&var"       as subg_var length=32,
            "&sgrow"     as sgrow length=8,
            cg.switch_seq,
            cg.SUBG,
            cats("&var","=",strip(put(cg.SUBG,best.))) as subg_group length=80,
            coalesce(r.n_units,0) as n_units
        from _cell_grid as cg
        left join _cell_raw as r
          on cg.switch_seq = r.switch_seq
         and cg.SUBG = r.SUBG
        order by cg.switch_seq, cg.SUBG;
    quit;

    proc sql noprint;
        select coalesce(min(n_units),0) into :minc from cell_cnt;
        select count(*) into :n_cell from cell_cnt;
    quit;

    data screen_one;
        length subg_var $32 data_src $12 sgrow $8 status $24;
        subg_var="&var"; data_src="&src_label"; sgrow="&sgrow";
        n_units=&n_unit; n_subjects=&n_subj; n_switch_seq=&n_seq; n_cell_rows=&n_cell; min_cell=&minc;
        if n_units=0 then status="NO_DATA";
        else if n_switch_seq ne 2 then status="NOT_2_SWITCH_SEQ";
        else if min_cell=0 then status="EMPTY";
        else if min_cell<10 then status="SPARSE";
        else status="OK";
    run;

    title "Subgroup count source: &var";
    proc print data=raw_freq noobs; run;
    proc print data=screen_one noobs; run;
    proc print data=cell_cnt noobs;
        var subg_var switch_seq SUBG subg_group n_units;
    run;
    title;

    %if %sysevalf(&n_unit>0 and &n_seq=2 and &minc>0) %then %do;
        ods exclude all;

        ods output FitStatistics=fit_cs Tests3=typ_cs Estimates=est_cs CovParms=cov_cs;
        proc mixed data=ana_long method=reml covtest noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER ACTIVE_SUBG CARRYOVER_SUBG/ solution ddfm=kr cl;
            repeated &tvar / subject=&unit type=cs;
            estimate "Trt SUBG0" ACTIVE 1 / cl;
            estimate "Trt SUBG1" ACTIVE 1 ACTIVE_SUBG 1 / cl;
            estimate "Trt int" ACTIVE_SUBG 1 / cl;
            estimate "Cry SUBG0" CARRYOVER 1 / cl;
            estimate "Cry SUBG1" CARRYOVER 1 CARRYOVER_SUBG 1 / cl;
            estimate "Cry int" CARRYOVER_SUBG 1 / cl;
        run;

        ods output FitStatistics=fit_un Tests3=typ_un Estimates=est_un CovParms=cov_un;
        proc mixed data=ana_long method=reml covtest noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER ACTIVE_SUBG CARRYOVER_SUBG/ solution ddfm=kr cl;
            repeated &tvar / subject=&unit type=un;
            estimate "Trt SUBG0" ACTIVE 1 / cl;
            estimate "Trt SUBG1" ACTIVE 1 ACTIVE_SUBG 1 / cl;
            estimate "Trt int" ACTIVE_SUBG 1 / cl;
            estimate "Cry SUBG0" CARRYOVER 1 / cl;
            estimate "Cry SUBG1" CARRYOVER 1 CARRYOVER_SUBG 1 / cl;
            estimate "Cry int" CARRYOVER_SUBG 1 / cl;
        run;

        ods output FitStatistics=fit_csn Tests3=typ_csn Estimates=est_csn CovParms=cov_csn;
        proc mixed data=ana_long method=reml covtest noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE CARRYOVER / solution ddfm=kr cl;
            repeated &tvar / subject=&unit type=cs;
            estimate "Reduced trt" ACTIVE 1 / cl;
            estimate "Reduced cry" CARRYOVER 1 / cl;
        run;

        ods output FitStatistics=fit_unn Tests3=typ_unn Estimates=est_unn CovParms=cov_unn;
        proc mixed data=ana_long method=reml covtest noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE CARRYOVER / solution ddfm=kr cl;
            repeated &tvar / subject=&unit type=un;
            estimate "Reduced trt" ACTIVE 1 / cl;
            estimate "Reduced cry" CARRYOVER 1 / cl;
        run;

        ods output FitStatistics=mlfit_cs_full;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER ACTIVE_SUBG CARRYOVER_SUBG / solution;
            repeated &tvar / subject=&unit type=cs;
        run;
        ods output FitStatistics=mlfit_cs_noact;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER CARRYOVER_SUBG / solution;
            repeated &tvar / subject=&unit type=cs;
        run;
        ods output FitStatistics=mlfit_cs_nocry;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER ACTIVE_SUBG / solution;
            repeated &tvar / subject=&unit type=cs;
        run;
        ods output FitStatistics=mlfit_cs_Reduced;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE CARRYOVER / solution;
            repeated &tvar / subject=&unit type=cs;
        run;

        ods output FitStatistics=mlfit_un_full;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER ACTIVE_SUBG CARRYOVER_SUBG / solution;
            repeated &tvar / subject=&unit type=un;
        run;
        ods output FitStatistics=mlfit_un_noact;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER CARRYOVER_SUBG / solution;
            repeated &tvar / subject=&unit type=un;
        run;
        ods output FitStatistics=mlfit_un_nocry;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE SUBG CARRYOVER ACTIVE_SUBG / solution;
            repeated &tvar / subject=&unit type=un;
        run;
        ods output FitStatistics=mlfit_un_Reduced;
        proc mixed data=ana_long method=ml noclprint;
            class &unit &tvar;
            model analysis_y=BASELINE PERIOD2 ACTIVE CARRYOVER / solution;
            repeated &tvar / subject=&unit type=un;
        run;
        ods exclude none;

        %tag_append(ds=fit_cs,out=fit_all,var=&var,model=CS,mean=full,src=&src_label);
        %tag_append(ds=fit_un,out=fit_all,var=&var,model=UN,mean=full,src=&src_label);
        %tag_append(ds=fit_csn,out=fit_all,var=&var,model=CSN,mean=Reduced,src=&src_label);
        %tag_append(ds=fit_unn,out=fit_all,var=&var,model=UNN,mean=Reduced,src=&src_label);
        %tag_append(ds=est_cs,out=est_all,var=&var,model=CS,mean=full,src=&src_label);
        %tag_append(ds=est_un,out=est_all,var=&var,model=UN,mean=full,src=&src_label);
        %tag_append(ds=est_csn,out=est_all,var=&var,model=CSN,mean=Reduced,src=&src_label);
        %tag_append(ds=est_unn,out=est_all,var=&var,model=UNN,mean=Reduced,src=&src_label);
        %tag_append(ds=typ_cs,out=typ_all,var=&var,model=CS,mean=full,src=&src_label);
        %tag_append(ds=typ_un,out=typ_all,var=&var,model=UN,mean=full,src=&src_label);
        %tag_append(ds=typ_csn,out=typ_all,var=&var,model=CSN,mean=Reduced,src=&src_label);
        %tag_append(ds=typ_unn,out=typ_all,var=&var,model=UNN,mean=Reduced,src=&src_label);
        %tag_append(ds=cov_cs,out=cov_all,var=&var,model=CS,mean=full,src=&src_label);
        %tag_append(ds=cov_un,out=cov_all,var=&var,model=UN,mean=full,src=&src_label);
        %tag_append(ds=cov_csn,out=cov_all,var=&var,model=CSN,mean=Reduced,src=&src_label);
        %tag_append(ds=cov_unn,out=cov_all,var=&var,model=UNN,mean=Reduced,src=&src_label);

        %get_ll(name=cs_full); %get_ll(name=cs_noact); %get_ll(name=cs_nocry); %get_ll(name=cs_Reduced);
        %get_ll(name=un_full); %get_ll(name=un_noact); %get_ll(name=un_nocry); %get_ll(name=un_Reduced);
        %make_lrt(var=&var,src=&src_label,cov=CS,test=Active interaction,full=cs_full,reduced=cs_noact,df=1);
        %make_lrt(var=&var,src=&src_label,cov=CS,test=Carryover interaction,full=cs_full,reduced=cs_nocry,df=1);
        %make_lrt(var=&var,src=&src_label,cov=CS,test=full versus Reduced,full=cs_full,reduced=cs_Reduced,df=3);
        %make_lrt(var=&var,src=&src_label,cov=UN,test=Active interaction,full=un_full,reduced=un_noact,df=1);
        %make_lrt(var=&var,src=&src_label,cov=UN,test=Carryover interaction,full=un_full,reduced=un_nocry,df=1);
        %make_lrt(var=&var,src=&src_label,cov=UN,test=full versus Reduced,full=un_full,reduced=un_Reduced,df=3);

        %if %sysfunc(exist(est_all)) %then %do;
            data key_est;
                set est_all;
                if mean_model="full" and Label in
                   ("Trt SUBG0","Trt SUBG1","Trt int","Cry SUBG0","Cry SUBG1","Cry int")
                then output;
                if mean_model="Reduced" and Label in ("Reduced trt","Reduced cry")
                then output;
            run;
            proc sort data=key_est; by mean_model model Label; run;
        %end;

        %if %sysfunc(exist(typ_all)) %then %do;
            data key_typ;
                set typ_all;
                where (mean_model="full" and Effect in ("ACTIVE_SUBG","CARRYOVER_SUBG"))
                   or (mean_model="Reduced" and Effect in ("ACTIVE","CARRYOVER"));
            run;
            proc sort data=key_typ; by mean_model model Effect; run;
        %end;

        title "Estimates: &var";
        %if %sysfunc(exist(key_est)) %then %do;
            proc print data=key_est noobs;
                var subg_var mean_model model Label Estimate StdErr DF tValue Probt Lower Upper;
            run;
        %end;

        title "Type III: &var";
        %if %sysfunc(exist(key_typ)) %then %do;
            proc print data=key_typ noobs;
                var subg_var mean_model model Effect NumDF DenDF FValue ProbF;
            run;
        %end;

        title "LRT: &var";
        %if %sysfunc(exist(lrt_all)) %then %do;
            proc print data=lrt_all noobs;
                var subg_var Covariance Test neg2ll_reduced neg2ll_full LR_Stat DF P_Value;
            run;
        %end;

        title "Fit stats: &var";
        %if %sysfunc(exist(fit_all)) %then %do;
            proc print data=fit_all noobs; run;
        %end;

        title;
    %end;
    %else %do;
        %put NOTE: Skip &var due to sparse cells. n_unit=&n_unit n_seq=&n_seq min_cell=&minc;
    %end;
%mend;


data var_map;
    length var $32 pos $20 neg $40;
    infile datalines dlm='|' truncover;
    input var :$32. pos :$20. neg :$40.;
    datalines;
TOBAC100|1|0
;
run;

data _null_;
    set var_map;
    call execute(cats('%nrstr(%fit_mcombm)(var=',strip(var),
                      ', pos=',strip(pos),
                      ', neg=',strip(neg),
                      ', sgrow=BASE);'));
run;

title;


/*============================================================*
 | Baseline outcome boxplots for 4 cells:
 | switch group x subgroup
 | Outcome: baseline CDRSUM = analysis_y
 *============================================================*/

%let subgvar = TOBAC100;
%let pos = 1;
%let neg = 0;

/* 1. Baseline-only dataset with subgroup coding */
data baseline_4cell;
    set cohort;
    where BASELINE = 1 and not missing(analysis_y);

    length SUBG_LABEL $30 CELL $80;

    if &subgvar in (&neg) then do;
        SUBG = 0;
        SUBG_LABEL = cats("&subgvar", "=0");
    end;
    else if &subgvar in (&pos) then do;
        SUBG = 1;
        SUBG_LABEL = cats("&subgvar", "=1");
    end;
    else delete;

    CELL = cats(switch_seq, " / ", SUBG_LABEL);

    keep NACCID tripleid switch_seq direction SUBG SUBG_LABEL CELL analysis_y;
run;

/* 2. QC: sample size and baseline outcome summary in each of 4 cells */
title "Baseline CDRSUM Summary by Switch Sequence and Subgroup";
proc sql;
    select 
        switch_seq,
        SUBG,
        SUBG_LABEL,
        count(*) as n,
        mean(analysis_y) as mean_CDRSUM format=8.3,
        std(analysis_y) as sd_CDRSUM format=8.3,
        median(analysis_y) as median_CDRSUM format=8.3,
        min(analysis_y) as min_CDRSUM format=8.3,
        max(analysis_y) as max_CDRSUM format=8.3
    from baseline_4cell
    group by switch_seq, SUBG, SUBG_LABEL
    order by switch_seq, SUBG;
quit;
title;


/*============================================================*
 | Figure
 | Baseline CDR-SB by Sequence and Subgroup
 |
 | Above each box:
 |   n = XX
 |   Mean (SD) = X.XX (X.XX)
 |
 | Below each box:
 |   Sequence, subgroup
 |   Range: Min-Max
 |
 | Boxplot:
 |   Diamond = Mean
 |   Horizontal line inside box = Median
 *============================================================*/


/*============================================================*
 | 1. Prepare plotting data
 *============================================================*/

data baseline_4cell_plot;
    set baseline_4cell;

    length seq_label  $12
           subg_label $8
           cell_label $30;

    /*--------------------------------------------------------*
     | Rename sequence
     *--------------------------------------------------------*/

    if switch_seq = "M_to_C" then
        seq_label = "AB sequence";

    else if switch_seq = "C_to_M" then
        seq_label = "BA sequence";

    else
        seq_label = strip(switch_seq);


    /*--------------------------------------------------------*
     | Rename subgroup
     *--------------------------------------------------------*/

    if SUBG = 0 then
        subg_label = "G=0";

    else if SUBG = 1 then
        subg_label = "G=1";

    else
        subg_label =
            cats(
                "G=",
                strip(
                    put(SUBG, best.)
                )
            );


    /*--------------------------------------------------------*
     | Manuscript cell label
     *--------------------------------------------------------*/

    cell_label =
        cats(
            seq_label,
            ", ",
            subg_label
        );


    /*--------------------------------------------------------*
     | Numeric x-axis position
     |
     | 1 = AB, G=0
     | 2 = AB, G=1
     | 3 = BA, G=0
     | 4 = BA, G=1
     *--------------------------------------------------------*/

    if seq_label = "AB sequence" and SUBG = 0 then
        cell_order = 1;

    else if seq_label = "AB sequence" and SUBG = 1 then
        cell_order = 2;

    else if seq_label = "BA sequence" and SUBG = 0 then
        cell_order = 3;

    else if seq_label = "BA sequence" and SUBG = 1 then
        cell_order = 4;

run;


/*============================================================*
 | 2. Calculate summary statistics
 *============================================================*/

proc summary data=baseline_4cell_plot nway;

    class
        cell_order
        cell_label;

    var
        analysis_y;

    output out=box_stats(
        drop=_TYPE_ _FREQ_
    )
        n      = n
        mean   = mean_CDRSB
        std    = sd_CDRSB
        median = median_CDRSB
        min    = min_CDRSB
        max    = max_CDRSB;

run;


/*============================================================*
 | 3. Optional QC table
 *============================================================*/

title "Baseline CDR-SB Summary by Sequence and Subgroup";

proc print data=box_stats noobs;

    var
        cell_order
        cell_label
        n
        mean_CDRSB
        sd_CDRSB
        median_CDRSB
        min_CDRSB
        max_CDRSB;

    format
        mean_CDRSB   8.3
        sd_CDRSB     8.3
        median_CDRSB 8.3
        min_CDRSB    8.3
        max_CDRSB    8.3;

run;

title;


/*============================================================*
 | 4. Create two-line x-axis labels
 |
 | Example:
 |
 | AB sequence, G=0
 | Range: 0.50-9.00
 |
 | "|" will later be used as line-break character.
 *============================================================*/

data box_stats_plot;
    set box_stats;

    length axis_label $100;

    axis_label =
        cats(
            strip(cell_label),

            "|Range: ",

            strip(
                put(min_CDRSB, 5.2)
            ),

            "-",

            strip(
                put(max_CDRSB, 5.2)
            )
        );

run;


/*============================================================*
 | 5. Create numeric format for x-axis
 |
 | cell_order=1 -> two-line AB/G=0 label
 | cell_order=2 -> two-line AB/G=1 label
 | etc.
 *============================================================*/

data cellfmt;
    set box_stats_plot;

    length
        fmtname $8
        type    $1
        label   $100;

    fmtname = "cellfmt";
    type    = "N";

    start = cell_order;
    end   = cell_order;

    label = axis_label;

    keep
        fmtname
        type
        start
        end
        label;

run;


proc format cntlin=cellfmt;
run;


/*============================================================*
 | 6. SG Annotation data set
 |
 | IMPORTANT:
 |
 | X position:
 |   X1SPACE = DATAVALUE
 |   X1 = 1, 2, 3, 4
 |
 | Y position:
 |   Y1SPACE = WALLPERCENT
 |
 | Therefore annotation stays INSIDE the plotting wall
 | and cannot move into the graph title.
 *============================================================*/

data anno_stats;
    set box_stats;

    length
        function   $9
        x1space    $12
        y1space    $12
        label      $100
        anchor     $12
        textcolor  $12
        textfont   $40
        textweight $12
        widthunit  $12;


    /*--------------------------------------------------------*
     | Common annotation properties
     *--------------------------------------------------------*/

    function = "text";

    /* X follows data category positions */
    x1space = "datavalue";

    /* Y follows percent of plot-wall height */
    y1space = "wallpercent";

    x1 = cell_order;

    anchor = "center";

    textcolor = "black";

    textfont = "Times New Roman";

    textweight = "normal";

    textsize = 8.5;


    /*--------------------------------------------------------*
     | Prevent long Mean (SD) text from wrapping
     *--------------------------------------------------------*/

    width = 24;

    widthunit = "percent";


    /*========================================================*
     | Annotation line 1:
     | Sample size
     |
     | 94% of plotting-wall height
     *========================================================*/

    y1 = 95;

    label =
        cats(
            "n = ",
            strip(
                put(n, 3.)
            )
        );

    output;


    /*========================================================*
     | Annotation line 2:
     | Mean (SD)
     |
     | 87% of plotting-wall height
     *========================================================*/

    y1 = 90;

    label =
        cats(
            "Mean (SD) = ",
            strip(
                put(mean_CDRSB, 5.2)
            ),
            " (",
            strip(
                put(sd_CDRSB, 5.2)
            ),
            ")"
        );

    output;


    keep
        function
        x1space
        y1space
        x1
        y1
        label
        anchor
        textcolor
        textfont
        textweight
        textsize
        width
        widthunit;

run;


/*============================================================*
 | 7. Optional QC:
 | Should contain exactly 8 rows:
 |
 | 4 groups x 2 annotation lines
 *============================================================*/

title "QC: SG Annotation Data";

proc print data=anno_stats noobs;
run;

title;


/*============================================================*
 | 8. Final plotting data
 |
 | Apply the custom two-line x-axis format.
 *============================================================*/

data baseline_4cell_final;
    set baseline_4cell_plot;

    format
        cell_order cellfmt.;

run;


/*============================================================*
 | 9. Sort in manuscript order
 *============================================================*/

proc sort data=baseline_4cell_final;
    by cell_order;
run;


/*============================================================*
 | 10. High-resolution output settings
 *============================================================*/

ods listing
    gpath="/home/u63719145/Clinical Paper/SubgroupEffect/real data"
   /*image_dpi=600 */
    ;


ods graphics on /
    reset

    width=7in

    /*
       Slightly taller than original figure to accommodate:
       - top statistics
       - two-line x-axis labels
    */
    height=6.2in

    imagename="Figure4(2)"

    outputfmt=svg ; /*svg OR png*/




/*============================================================*
 | 12. Final boxplot
 *============================================================*/

proc sgplot
    data=baseline_4cell_final
    sganno=anno_stats
    noautolegend;


    /*========================================================*
     | Four boxplots
     *========================================================*/

    vbox analysis_y /
        category=cell_order;


    /*========================================================*
     | X axis
     |
     | Example:
     |
     | AB sequence, G=0
     | Range: 0.50-9.00
     *========================================================*/

    xaxis

        label=" "

        values=(
            1
            2
            3
            4
        )

        fitpolicy=splitalways

        splitchar="|"

        labelattrs=(
            family="Times New Roman"
            size=10pt
        )

        valueattrs=(
            family="Times New Roman"
            size=8.5pt
        );


    /*========================================================*
     | Y axis
     |
     | Actual maximum = 18.
     |
     | max=24 intentionally leaves blank space above
     | the observations.
     |
     |
     | The 18-point outlier is approximately 75% of wall
     | height, leaving separation from annotations.
     *========================================================*/

    yaxis

        label="Baseline CDR-SB"

        min=0

        max=28

        values=(
            0
            5
            10
            15
            20
            25
        )

        labelattrs=(
            family="Times New Roman"
            size=10pt
        )

        valueattrs=(
            family="Times New Roman"
            size=9pt
        );

run;


/*============================================================*
 | 13. Finish output
 *============================================================*/

title;

ods graphics off;

ods listing close;


