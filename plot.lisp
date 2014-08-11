(in-package :cl-nsga2)

(defun plot-front (pop &key title filename (stream t) (if-exists :supersede)
			 (obj1 0) (obj2 1))
  "Plot Pareto front"
  (flet ((write-plot (s)
	   (when title (format s "set title \"~a\"~%" title))
	   (format s "plot '-' with points~&")
	   (dolist (ind pop)
	     (with-slots (obj) ind
	       (format s "~f ~f~%" (elt obj obj1) (elt obj obj2))))))
    (if filename ; write to file if specified
	(with-open-file (s filename :direction :output :if-exists if-exists)
	  (write-plot s))
	(write-plot stream)))) ; otherwise write to stream
