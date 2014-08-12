(ql:quickload :cl-nsga2)
(ql:quickload :bld-ode)
(ql:quickload :bld-utils)

(in-package :cl-nsga2)

(use-package :bld-ode)
(use-package :bld-utils)

(defparameter *sail-circ-coplanar-param*
  (make-hash*
   :be 0.1d0
   :mu 1d0
   :x0 (list 1d0 0d0 0d0 1d0)
   :xf (list 1.5d0 0d0 0d0 (sqrt (/ mu 1.5d0)))
   :t0 0d0
   :tf 10d0))

(defun sail-circ-coplanar-eom (tm x u)
  (lethash (be mu) *sail-circ-coplanar-param*
    (destructuring-bind (r th vr vt) x
      (let ((co (cos u))
	    (si (sin u)))
	(list vr
	      (/ vt r)
	      (+ (/ (expt vt 2) r)
		 (* mu (/ (1- (* be (expt co 2) (abs co))) (expt r 2))))
	      (- (/ (* mu be (expt co 2) si) (expt r 2))
		 (/ (* vr vt) r)))))))

(lethash (x0 t0) *sail-circ-coplanar-param*
  (defun sail-circ-coplanar-prop (tf utm)
    (let* ((nvar (length utm))
	   (dtf (/ (- tf t0) nvar)))
      (loop for u in utm
	 for t0i = t0 then tfi
	 for tfi = (+ t0i dtf)
	 for x0i = x0 then xfi
	 for traji = (rka #'sail-circ-coplanar-eom t0i tfi x0i :hmax (- tfi t0i) :param u)
	 for xfi = (second (car (last traji)))
	 append traji))))

(lethash (xf tf) *sail-circ-coplanar-param*
  (defun sail-circ-coplanar-objfun (ind)
    (let* ((trajectory (sail-circ-coplanar-prop tf (xvar ind)))
	   (xftraj (second (car (last trajectory))))
	   (xferr (list (abs (- (elt xftraj 0) (elt xf 0)))
			(abs (- (elt xftraj 2) (elt xf 2)))
			(abs (- (elt xftraj 3) (elt xf 3))))))
      (values xferr nil))))

(defparameter *sail-circ-coplanar-options*
  (make-instance
   'options
   :popsize 200
   :ngen 100
   :nobj 3
   :ncon 0
   :nvar 10
   :minvar (make-list 10 :initial-element (- (/ pi 2d0)))
   :maxvar (make-list 10 :initial-element (/ pi 2d0))
   :pcross 0.9d0
   :pmut 0.5d0
   :eta-c 10d0
   :eta-m 20d0
   :objfun #'sail-circ-coplanar-objfun))

(defun sail-circ-coplanar-plot (ind &key filename title (stream t) (if-exists :supersede))
  (flet ((sail-plot (s)
	   (when title (format s "set title \"~a\"~%" title))
	   (format s "set polar~%")
	   (format s "set size ratio -1~%")
	   (format s "set grid polar~%")
	   (format s "plot '-' with lines title \"\"~%")
	   (loop with tf = (car (xvar ind))
	      with utm = (cdr (xvar ind))
	      for (tm x) in (sail-circ-coplanar-prop tf utm)
	      for (r th vr vt) = x
	      do (format s "~f ~f~%" th r))))
    (if filename
	(with-open-file (s filename :direction :output :if-exists if-exists)
	  (sail-plot s))
	(sail-plot stream))))

;;; Optimize for time of flight

(lethash (xf) *sail-circ-coplanar-tf-param*
  (defun sail-circ-coplanar-tf-objfun (ind)
    (with-slots (xvar) ind
      (let* ((tf (car xvar))
	     (utm (cdr xvar))
	     (trajectory (sail-circ-coplanar-prop tf utm)) ; prop
	     (xftraj (second (car (last trajectory)))) ; XF of traj
	     (xferr (list (abs (- (elt xftraj 0) (elt xf 0)))
			  (abs (- (elt xftraj 2) (elt xf 2)))
			  (abs (- (elt xftraj 3) (elt xf 3))))))
	(values
	 (list tf ; TF
	       (norminfx xferr))
	 nil)))))

(defparameter *sail-circ-coplanar-tf-options*
  (make-instance
   'options
   :popsize 200
   :ngen 400
   :nobj 2
   :ncon 0
   :nvar 11
   :minvar (cons 5 (make-list 10 :initial-element (- (/ pi 2d0))))
   :maxvar (cons 20 (make-list 10 :initial-element (/ pi 2d0)))
   :pcross 0.9d0
   :pmut 0.5d0
   :eta-c 10d0
   :eta-m 20d0
   :objfun #'sail-circ-coplanar-tf-objfun))

(defun run-sail-tf-problem ()
  (defparameter *sail-circ-coplanar-tf-pop*
    (nsga2 *sail-circ-coplanar-tf-options*))
  (sail-circ-coplanar-plot 
   (find-min *sail-circ-coplanar-tf-pop* #'(lambda (ind) (cdr (obj ind))))
   :filename "plots/sail-circ-coplanar-tf.gp"
   :title "Circular coplanar minimum time sail transfer from Earth to Mars")
  (plot-front
   *sail-circ-coplanar-tf-pop*
   :title "Pareto front of minimum time circular coplanar sail transfer from Earth to Mars"
   :filename "plots/sail-circ-coplanar-tf-front.gp"))
