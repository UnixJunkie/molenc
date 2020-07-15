
module Fn = Filename

let plot_histograms histo_data_fn mean1 mean2 =
  let tmp_gnuplot_fn = Fn.temp_file "gnuplot_" ".gpl" in
  Utls.string_list_to_file tmp_gnuplot_fn
    ["set xlabel 'score'";
     "set ylabel 'frequency'";
     "gauss1(x) = a1/(sqrt(2*pi)*sigma1)*exp(-(x-mean1)**2/(2*sigma1**2))";
     "gauss2(x) = a2/(sqrt(2*pi)*sigma2)*exp(-(x-mean2)**2/(2*sigma2**2))";
     sprintf "mean1 = %d" mean1;
   sprintf "mean2 = %d" mean2;
     sprintf "fit gauss1(x) '%s' u 1:2 via a1,sigma1,mean1" histo_data_fn;
     sprintf "fit gauss2(x) '%s' u 1:3 via a2,sigma2,mean2" histo_data_fn;
     sprintf "plot '%s' u 1:2 w l t 'smaller sample', \\" histo_data_fn;
     "''        u 1:3 w l t 'bigger sample', \\";
     "gauss1(x) t 'smaller fit', \\";
     "gauss2(x) t 'bigger fit'"];
  Log.info "Gnuplot: script: %s" tmp_gnuplot_fn;
  Utls.run_command (sprintf "gnuplot -persist %s" tmp_gnuplot_fn)
