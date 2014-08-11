;;; Test optimization problems

(in-package :cl-nsga2)

(defun ctp1-objfun (ind)
  (with-slots (xvar) ind
    (let* ((g (1+ (elt xvar 1)))
	   (obj0 (elt xvar 0))
	   (obj1 (* g (exp (- (/ obj0 g)))))
	   (constr0 (constraint-gte (/ obj1 (* 0.858d0 (exp (* -0.541d0 obj0)))) 1d0))
	   (constr1 (constraint-gte (/ obj1 (* 0.728d0 (exp (* -0.295d0 obj0)))) 1d0)))
      (values
       (list obj0 obj1)
       (list constr0 constr1)))))

(defparameter *ctp1-options*
  (make-instance 
   'options
   :popsize 200
   :ngen 100
   :nobj 2
   :ncon 2
   :nvar 2
   :minvar (list 0d0 0d0)
   :maxvar (list 1d0 1d0)
   :pcross 0.9d0
   :pmut 0.5d0
   :eta-c 10d0
   :eta-m 20d0
   :objfun #'ctp1-objfun))

(defun binh-korn-objfun (ind)
  (with-slots (xvar) ind
    (destructuring-bind (x y) xvar
      (let ((f1 (+ (* 4 (expt x 2))
		   (* 4 (expt y 2))))
	    (f2 (+ (expt (- x 5) 2)
		   (expt (- y 5) 2)))
	    (g1 (constraint-lte
		 (+ (expt (- x 5) 2) (expt y 2)) 25))
	    (g2 (constraint-gte
		 (+ (expt (- x 8) 2) (expt (+ y 3) 2)) 7.7)))
	(values (list f1 f2) (list g1 g2))))))

(defparameter *binh-korn-options*
  (make-instance 
   'options
   :popsize 200
   :ngen 100
   :nobj 2
   :ncon 2
   :nvar 2
   :minvar (list 0d0 0d0)
   :maxvar (list 5d0 3d0)
   :pcross 0.9d0
   :pmut 0.5d0
   :eta-c 10d0
   :eta-m 20d0
   :objfun #'binh-korn-objfun))

(defun fonseca-fleming-objfun (ind)
  (with-slots (xvar) ind
    (let ((n (length xvar)))
      (values
       (list (- 1 (exp (- (reduce #'+ (mapcar #'(lambda (xi) (expt (- xi (/ (sqrt n))) 2)) xvar)))))
	     (- 1 (exp (- (reduce #'+ (mapcar #'(lambda (xi) (expt (+ xi (/ (sqrt n))) 2)) xvar))))))
       nil))))

(defparameter *fonseca-fleming-options*
  (make-instance 
   'options
   :popsize 200
   :ngen 100
   :nobj 2
   :ncon 0
   :nvar 10
   :minvar (make-list 10 :initial-element -4d0)
   :maxvar (make-list 10 :initial-element 4d0)
   :pcross 0.9d0
   :pmut 0.5d0
   :eta-c 10d0
   :eta-m 20d0
   :objfun #'fonseca-fleming-objfun))

(ql:quickload :bld-ode)
(ql:quickload :bld-utils)
(use-package :bld-utils)

(defparameter *sail-circ-coplanar-param*
  (make-hash*
   :be 0.1d0
   :mu 1d0
   :x0 (list 1d0 0d0 0d0 1d0)
   :xf (list 1.5d0 0d0 0d0 (sqrt (/ mu 1.5d0)))
   :xftol (list 1d-3 1d-3 1d-3) ; theta not used
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

(lethash (x0 xf t0 tf) *sail-circ-coplanar-param*
  (defun sail-circ-coplanar-prop (xvar)
    (let* ((nvar (length xvar))
	   (dtf (/ (- tf t0) nvar)))
      (loop for u in xvar
	 for t0i = t0 then tfi
	 for tfi = (+ t0i dtf)
	 for x0i = x0 then xfi
	 for traji = (rka #'sail-circ-coplanar-eom t0i tfi x0i :hmax (- tfi t0i) :param u)
	 for xfi = (second (car (last traji)))
	 append traji))))

(lethash (xf xftol) *sail-circ-coplanar-param*
  (defun sail-circ-coplanar-objfun (ind)
    (let* ((trajectory (sail-circ-coplanar-prop (xvar ind)))
	   (xftraj (second (car (last trajectory))))
	   (objvals
	    (list
	     (abs (- (first xftraj) (first xf)))
	     (abs (- (third xftraj) (third xf)))
	     (abs (- (fourth xftraj) (fourth xf)))))
	   (constr
	    (list
	     (abs (- (first objvals) (first xftol)))
	     (abs (- (second objvals) (second xftol)))
	     (abs (- (third objvals) (third xftol))))))
      (values objvals constr))))

(defparameter *sail-circ-coplanar-options*
  (make-instance
   'options
   :popsize 200
   :ngen 100
   :nobj 3
   :ncon 3
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
	   (loop for (tm x) in (sail-circ-coplanar-prop
				(xvar ind))
	      for (r th vr vt) = x
	      do (format s "~f ~f~%" th r))))
    (if filename
	(with-open-file (s filename :direction :output :if-exists if-exists)
	  (sail-plot s))
	(sail-plot stream))))

