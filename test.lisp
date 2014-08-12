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

(defun chakong-haimes-objfun (ind)
  (with-slots (xvar) ind
    (destructuring-bind (x y) xvar
      (values
       (list
	(+ 2 (expt (- x 2) 2) (expt (- y 1) 2))
	(+ (* 9 x) (expt (- y 1) 2)))
       (list
	(constraint-lte (+ (expt x 2) (expt y 2)) 225)
	(constraint-lte (+ x (* -3 y) 10) 0))))))

(defparameter *chakong-haimes-options*
  (make-instance
   'options
   :popsize 200
   :ngen 100
   :nobj 2
   :ncon 2
   :nvar 2
   :minvar '(-20 -250)
   :maxvar '(250 20)
   :pcross 0.9d0
   :pmut 0.5d0
   :eta-c 10d0
   :eta-m 20d0
   :objfun #'chakong-haimes-objfun))

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

