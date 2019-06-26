%include 'All_macros.sas';

%let nboot = 1;
%let nboot_sens = 1;

/* Table 1 */
%Table1(nboot = &nboot, row = 1, eof = &eof);
%Table1(nboot = &nboot, row = 2, eof = &eof);

%Table1(nboot = &nboot, row = 3, eof = &eof);
%Table1(nboot = &nboot, row = 4, eof = &eof);
%Table1(nboot = &nboot, row = 5, eof = &eof);

/* Table 2*/
%Table2(nboot = &nboot, eof = &eof, row = 1, sens = 0);
%Table2(nboot = &nboot, eof = &eof, row = 2, sens = 0);
%Table2(nboot = &nboot, eof = &eof, row = 3, sens = 0);
%Table2(nboot = &nboot, eof = &eof, row = 4, sens = 0);
%Table2(nboot = &nboot, eof = &eof, row = 5, sens = 0);
%Table2(nboot = &nboot, eof = &eof, row = 6, sens = 0);
%Table2(nboot = &nboot, eof = &eof, row = 7, sens = 0);
%Table2(nboot = &nboot, eof = &eof, row = 8, sens = 0);

/* Table 2 sensitivity analyses*/
%Table2(nboot = &nboot_sens, eof = &eof, row = 1, sens = 1, sens_option = "Model type");
%Table2(nboot = &nboot_sens, eof = &eof, row = 2, sens = 1, sens_option = "Model type");
%Table2(nboot = &nboot_sens, eof = &eof, row = 1, sens = 1, sens_option = "Adherence missingness");
%Table2(nboot = &nboot_sens, eof = &eof, row = 4, sens = 1, sens_option = "Censoring time");
%Table2(nboot = &nboot_sens, eof = &eof, row = 5, sens = 1, sens_option = "Censoring time");
%Table2(nboot = &nboot_sens, eof = &eof, row = 6, sens = 1, sens_option = "Censoring time");
%Table2(nboot = &nboot_sens, eof = &eof, row = 7, sens = 1, sens_option = "Censoring time");
%Table2(nboot = &nboot_sens, eof = &eof, row = 8, sens = 1, sens_option = "Censoring time");
%Table2(nboot = &nboot_sens, eof = &eof, row = 6, sens = 1, sens_option = "IP weights");
%Table2(nboot = &nboot_sens, eof = &eof, row = 8, sens = 1, sens_option = "IP weights");
%Table2(nboot = &nboot_sens, eof = &eof, row = 7, sens = 1, sens_option = "Functional form");
%Table2(nboot = &nboot_sens, eof = &eof, row = 8, sens = 1, sens_option = "Functional form");

/* Table 3*/
%Table3(nboot = &nboot, eof = &eof, effect = "ITT", adjust = "covs0", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "none", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "substudy", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "covs0", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "covs_t", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "dose-response", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "none", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "substudy", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "covs0", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "covs_t", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "dose-response", est = "HR", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "none", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "substudy", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "covs0", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "covs_t", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "dose-response", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "none", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "substudy", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "covs0", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "covs_t", est = "RD", sens = 0);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "dose-response", est = "RD", sens = 0);


/* Table 3 sensitivity analyses*/
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, major AE", adjust = "covs_t", est = "HR", sens = 1, sens_option = "Censoring time");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, major AE", adjust = "dose-response", est = "HR", sens = 1, sens_option = "Censoring time");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, major AE", adjust = "covs_t", est = "HR", sens = 1, sens_option = "IP weights");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, major AE", adjust = "dose-response", est = "HR", sens = 1, sens_option = "IP weights");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, major AE", adjust = "dose-response", est = "HR", sens = 1, sens_option = "Functional form");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, any AE", adjust = "covs_t", est = "HR", sens = 1, sens_option = "Censoring time");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, any AE", adjust = "dose-response", est = "HR", sens = 1, sens_option = "Censoring time");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, any AE", adjust = "covs_t", est = "HR", sens = 1, sens_option = "IP weights");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, any AE", adjust = "dose-response", est = "HR", sens = 1, sens_option = "IP weights");
%Table3(nboot = &nboot_sens, eof = &eof, effect = "PPE, any AE", adjust = "dose-response", est = "HR", sens = 1, sens_option = "Functional form");

/* Figure 1*/
%Table3(nboot = &nboot, eof = &eof, effect = "ITT", adjust = "covs0", est = "RD", sens = 0, graph = 1);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, major AE", adjust = "covs_t", est = "RD", sens = 0, graph = 1);
%Table3(nboot = &nboot, eof = &eof, effect = "PPE, any AE", adjust = "covs_t", est = "RD", sens = 0, graph = 1);
