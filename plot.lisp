(in-package :cl-nsga2)

(defun plot-front (pop &key title filename (stream t) (if-exists :supersede)
			 (obj1f #'first) (obj2f #'(lambda (obj) (normelist (cdr obj)))))
  "Plot Pareto front"
  (flet ((write-plot (s)
	   (when title (format s "set title \"~a\"~%" title))
	   (format s "plot '-' with points title \"\"~&")
	   (dolist (ind pop)
	     (with-slots (obj) ind
	       (format s "~f ~f~%" (funcall obj1f obj) (funcall obj2f obj))))))
    (if filename ; write to file if specified
	(with-open-file (s filename :direction :output :if-exists if-exists)
	  (write-plot s))
	(write-plot stream)))) ; otherwise write to stream
