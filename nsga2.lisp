(in-package :cl-nsga2)

(defparameter *inf* 1d14)

(defparameter *eps* 1d-14)

(defclass options ()
  ((popsize :initarg :popsize)
   (ngen :initarg :ngen)
   (nobj :initarg :nobj)
   (ncon :initarg :ncon)
   (nvar :initarg :nvar)
   (minvar :initarg :minvar)
   (maxvar :initarg :maxvar)
   (pcross :initarg :pcross)
   (pmut :initarg :pmut)
   (eta-c :initarg :eta-c)
   (eta-m :initarg :eta-m)
   (objfun :initarg :objfun)))

(defclass individual ()
  ((rank :accessor rank)
   (constr-sum :accessor constr-sum)
   (constr-count :accessor constr-count)
   (xvar :initarg :xvar :accessor xvar)
   (obj :initarg :obj :accessor obj)
   (constr :accessor constr)
   (crowd-dist :accessor crowd-dist)
   (s-dom :accessor s-dom)
   (n-dom :accessor n-dom)))

(defun random-range (low high)
  (+ low (random (- high low))))

(defun initialize-ind (options)
  (with-slots (minvar maxvar) options
    (make-instance 'individual
     :xvar (loop for minv in minvar
	      for maxv in maxvar
	      collect (random-range minv maxv)))))

(defun initialize-pop (options)
  (with-slots (popsize) options
    (loop repeat popsize
       collect (initialize-ind options))))

(defun evaluate-ind (options ind)
  (with-slots (objfun ncon) options
    (with-slots (obj constr constr-sum constr-count) ind
      (multiple-value-bind (obj-f constr-f) (funcall objfun ind)
	(setf obj obj-f constr constr-f)
	(setf constr-sum 0d0)
	(setf constr-count 0)
	(loop for constrj in constr
	   when (< constrj 0)
	   do (incf constr-sum constrj)
	     (incf constr-count)))))
  ind)

(defun evaluate-pop (options pop)
  (mapcar #'(lambda (ind) (evaluate-ind options ind)) pop))

(defun dominates (p q)
  (with-slots ((sum-p constr-sum) (count-p constr-count) (obj-p obj)) p
    (with-slots ((sum-q constr-sum) (count-q constr-count) (obj-q obj)) q
      (cond
	;; p and q feasible
	((and (zerop count-p) (zerop count-q))
	 (let (pdomq qdomp)
	   (loop for op in obj-p
	      for oq in obj-q
	      if (< op oq)
	      do (setf pdomq t)
	      else if (> op oq)
	      do (setf qdomp t))
	   (cond
	     ((and pdomq (not qdomp)) p)
	     ((and (not pdomq) qdomp) q))))
	;; p feasible, q infeasible
	((and (zerop count-p) (not (zerop count-q))) p)
	;; q feasible, p infeasible
	((and (not (zerop count-p)) (zerop count-q)) q)
	;; p and q infeasible
	(t (cond
	     ((< sum-p sum-q) p)
	     ((< sum-p sum-q) q)))))))

(defun fast-non-dominated-sort (pop)
  "Sort test population into Pareto fronts, updating the number-dominating (n-dom) and set-of-dominated (s-dom) slots of each"
  ;; Initialize first front
  (let (f0)
    (dolist (p pop)
      (setf (s-dom p) nil)
      (setf (n-dom p) 0)
      (dolist (q pop)
	(let ((dom (dominates p q)))
	  (cond
	    ((eq dom p) (push q (s-dom p)))
	    ((eq dom q) (incf (n-dom p))))))
      (when (zerop (n-dom p))
	(setf (rank p) 0)
	(push p f0)))
    ;; Subsequent fronts
    (loop for front = 0 then (incf front)
       for fi = f0 then
	 (loop for p in fi
	    nconc (loop for q in (s-dom p)
		     do (decf (n-dom q))
		     when (zerop (n-dom q))
		     do (setf (rank q) front)
		     and collect q))
       while fi
       collect fi)))

(defun assign-crowding-distance-front (front &aux (numobj (length (obj (first front)))))
  (dolist (ind front) (setf (crowd-dist ind) 0d0))
  (loop for m below numobj
     do (setf front (sort front #'< :key #'(lambda (p) (nth m (obj p))))
	      (crowd-dist (first front)) *inf*
	      (crowd-dist (car (last front))) *inf*)
       (loop for i from 1 below (1- (length front))
	  do (incf (crowd-dist (elt front i))
		   (- (elt (obj (elt front (1+ i))) m)
		      (elt (obj (elt front (1- i))) m)))))
  front)

(defun assign-crowding-distance (fronts)
  (loop for front in fronts
     collect (assign-crowding-distance-front front)))

(defun crowd-compare (p q)
  (or (< (rank p) (rank q))
      (and (= (rank p) (rank q))
	   (> (crowd-dist p) (crowd-dist q)))))

(defun tournament (p q)
  (let ((dom (dominates p q)))
    (cond
      ((eq dom p) p)
      ((eq dom q) q)
      ((> (crowd-dist p) (crowd-dist q)) p)
      ((> (crowd-dist q) (crowd-dist p)) q)
      (t (if (<= (random 1d0) 0.5d0) p q)))))

(defun crossover (options parent1 parent2)
  (with-slots (pcross nvar minvar maxvar eta-c) options
    (if (<= (random 1d0) pcross) ; do crossover?
	(loop with child1 = (make-instance 'individual :xvar (make-list nvar :initial-element 0d0) :obj (make-list nvar :initial-element 0d0))
	   with child2 = (make-instance 'individual :xvar (make-list nvar :initial-element 0d0) :obj (make-list nvar :initial-element 0d0))
	   for x1 in (xvar parent1)
	   for x2 in (xvar parent2)
	   for yl in minvar
	   for yu in maxvar
	   for i below nvar
	   do (if (<= (random 1d0) 0.5d0) ; Even chance of copying parent1/2 to child1/2
		  (if (> (abs (- x1 x2)) *eps*) ; difference larger than *EPS*
		      (destructuring-bind (y1 y2)
			  (if (< x1 x2) ; smaller into y1, larger into y2
			      (list x1 x2)
			      (list x2 x1))
			(let* ((rand (random 1d0)) ; uniform random # to SBX PDF
			       (beta (1+ (/ (* 2d0 (- y1 yl)) (- y2 y1))))
			       (alpha (- 2d0 (expt beta (- (1+ eta-c)))))
			       (betaq (if (<= rand (/ alpha))
					  (expt (* rand alpha) (/ (1+ eta-c)))
					  (expt (/ (- 2d0 (* rand alpha))) (/ (1+ eta-c)))))
			       (c1 (* 0.5d0 (+ y1 y2 (- (* betaq (- y2 y1)))))) ; child 1 (lower)
			       (beta (1+ (* 2d0 (/ (- yu y2) (- y2 y1)))))
			       (alpha (- 2d0 (expt beta (- (1+ eta-c)))))
			       (betaq (if (<= rand (/ alpha))
					  (expt (* rand alpha) (/ (1+ eta-c)))
					  (expt (/ (- 2d0 (* rand alpha))) (/ (1+ eta-c)))))
			       (c2 (* 0.5d0 (+ y1 y2 (* betaq (- y2 y1))))) ; child 2 (upper)
			       (c1 (max (min c1 yu) yl)) ; apply limits
			       (c2 (max (min c2 yu) yl))) ; apply limits
			  (if (<= (random 1d0) 0.5d0) ; even chance of child 1 or 2 getting var
			      (setf (elt (obj child1) i) c2
				    (elt (obj child2) i) c1)
			      (setf (elt (obj child1) i) c1
				    (elt (obj child2) i) c2))))
		      (setf (elt (obj child1) i) x1 ; if |x1-x2| < *EPS* just copy
			    (elt (obj child2) i) x2))
		  (setf (elt (obj child1) i) x1 ; even chance of just copying parent var
			(elt (obj child2) i) x2))
	   finally (return (list child1 child2)))
	;; Else copy parents design variables into new children
	(list (make-instance 'individual :xvar (copy-list (xvar parent1)))
		(make-instance 'individual :xvar (copy-list (xvar parent2)))))))

(defun crossover2 (options parent1 parent2)
  (with-slots (pcross nvar minvar maxvar eta-c) options
    (let ((child1 (make-instance 'individual :xvar (make-list (length (xvar parent1)))))
	  (child2 (make-instance 'individual :xvar (make-list (length (xvar parent2)))))
	  i 
	  rand
	  y1 y2 yl yu
	  c1 c2
	  alpha beta betaq)
      (if (<= (random 1d0) pcross)
	  (dotimes (i nvar)
	    (if (<= (random 1d0) 0.5d0)
		(if (> (abs (- (elt (xvar parent1) i) (elt (xvar parent2) i))) *eps*)
		    (progn
		      (if (< (elt (xvar parent1) i) (elt (xvar parent2) i))
			  (setf y1 (elt (xvar parent1) i)
				y2 (elt (xvar parent2) i))
			  (setf y1 (elt (xvar parent2) i)
				y2 (elt (xvar parent1) i)))
		      (setf yl (elt minvar i))
		      (setf yu (elt maxvar i))
		      (setf rand (random 1d0))
		      (setf beta (+ 1d0 (/ (* 2d0 (- y1 yl)) (- y2 y1))))
		      (setf alpha (- 2d0 (expt beta (- (+ eta-c 1d0)))))
		      (if (<= rand (/ 1d0 alpha))
			  (setf betaq (expt (* rand alpha) (/ 1d0 (+ eta-c 1d0))))
			  (setf betaq (expt (/ 1d0 (- 2d0 (* rand alpha))) (/ 1d0 (+ eta-c 1d0)))))
		      (setf c1 (* 0.5d0 (- (+ y1 y2) (* betaq (- y2 y1)))))
		      (setf beta (+ 1d0 (/ (* 2d0 (- yu y2)) (- y2 y1))))
		      (setf alpha (- 2d0 (expt beta (- (+ eta-c 1d0)))))
		      (if (<= rand (/ 1d0 alpha))
			  (setf betaq (expt (* rand alpha) (/ 1d0 (+ eta-c 1d0))))
			  (setf betaq (expt (/ 1d0 (- 2d0 (* rand alpha))) (/ 1d0 (+ eta-c 1d0)))))
		      (setf c2 (* 0.5d0 (+ (+ y1 y2) (* betaq (- y2 y1)))))
		      (if (< c1 yl) (setf c1 yl))
		      (if (< c2 yl) (setf c2 yl))
		      (if (> c1 yu) (setf c1 yu))
		      (if (> c2 yu) (setf c2 yu))
		      (if (<= (random 1d0) 0.5d0)
			  (setf (elt (xvar child1) i) c2
				(elt (xvar child2) i) c1)
			  (setf (elt (xvar child1) i) c1
				(elt (xvar child2) i) c2)))
		    (setf (elt (xvar child1) i) (elt (xvar parent1) i)
			  (elt (xvar child2) i) (elt (xvar parent2) i)))
		(setf (elt (xvar child1) i) (elt (xvar parent1) i)
		      (elt (xvar child2) i) (elt (xvar parent2) i))))
	  (dotimes (i nvar)
	    (setf (elt (xvar child1) i) (elt (xvar parent1) i)
		  (elt (xvar child2) i) (elt (xvar parent2) i))))
      (list child1 child2))))

(defun selection (options parents)
  (with-slots (popsize) options
    (destructuring-bind (a1 a2)
	(loop for i below popsize collect i into a1 collect i into a2
	   finally (return (list a1 a2)))
      (dotimes (i popsize)
	(rotatef (elt a1 i) (elt a1 (random-range i popsize)))
	(rotatef (elt a2 i) (elt a2 (random-range i popsize))))
      (loop for i = 0 then (incf i 4)
	 while (< i popsize)
	 for parent11 = (tournament (elt parents (elt a1 i)) (elt parents (elt a1 (+ i 1))))
	 for parent21 = (tournament (elt parents (elt a1 (+ i 2))) (elt parents (elt a1 (+ i 3))))
	 for (child11 child21) = (crossover2 options parent11 parent21)
	 for parent12 = (tournament (elt parents (elt a2 i)) (elt parents (elt a2 (+ i 1))))
	 for parent22 = (tournament (elt parents (elt a2 (+ i 2))) (elt parents (elt a2 (+ i 3))))
	 for (child12 child22) = (crossover2 options parent12 parent22)
	 collect child11
	 collect child21
	 collect child12
	 collect child22))))

(defun mutation-ind (options ind)
  (with-slots (nvar pmut eta-m minvar maxvar) options
    (with-slots (xvar) ind
      (dotimes (j nvar)
	(when (<= (random 1d0) pmut)
	  (let* ((y (elt xvar j))
		 (yl (elt minvar j))
		 (yu (elt maxvar j))
		 (delta1 (/ (- y yl) (- yu yl)))
		 (delta2 (/ (- yu y) (- yu yl)))
		 (rnd (random 1d0))
		 (mut-pow (/ (1+ eta-m)))
		 (xy (if (<= rnd 0.5d0)
			 (- 1d0 delta1)
			 (- 1d0 delta2)))
		 (val (if (<= rnd 0.5d0)
			  (+ (* 2d0 rnd)
			     (* (- 1d0 (* 2d0 rnd)) (1+ (expt xy (1+ eta-m)))))
			  (+ (* 2d0 (- 1d0 rnd))
			     (* 2d0 (- rnd 0.5d0) (expt xy (1+ eta-m))))))
		 (deltaq (if (<= rnd 0.5d0)
			     (- (expt val mut-pow) 1d0)
			     (- 1d0 (expt val mut-pow))))
		 (y (incf y (* deltaq (- yu yl))))
		 (y (min (max y yl) yu)))
	    (setf (elt xvar j) y))))))
  ind)

(defun mutation-pop (options pop)
  (loop for ind in pop
     collect (mutation-ind options ind)))

(defun build-new-parent (fronts n)
  (loop for front in fronts
     append front into new-parents
     while (< (length new-parents) n)
     finally (return (subseq new-parents 0 n))))

(defun nsga2 (options)
  (with-slots (popsize ngen nobj ncon nvar minvar maxvar pcross pmut eta-c eta-m) options
    (let* ((parents (evaluate-pop 
		     options (initialize-pop options)))
	   (pfronts (assign-crowding-distance
		     (fast-non-dominated-sort parents))))
      (format t "1 ")
      (dotimes (i (1- ngen))
	(format t "~a " (+ i 2))
	(let* ((children (evaluate-pop 
		      options
		      (mutation-pop
		       options (selection options parents))))
	       (mixed (append parents children))
	       (mixed-fronts
		(assign-crowding-distance
		 (assign-crowding-distance
		  (fast-non-dominated-sort mixed))))
	       (mixed (sort mixed #'crowd-compare)))
	  (setf parents (subseq mixed 0 popsize))))
      (values parents))))

(defun norme2obj (ind)
  "Euclidean norm squared of objective function values"
  (reduce #'+ (mapcar #'(lambda (obj) (expt obj 2)) (obj ind))))

(defun normeobj (ind)
  "Euclidean norm of objective function values"
  (sqrt (norme2obj ind)))

(defun find-min-obj (pop)
  (loop for ind in pop
     for objnorm = (normeobj ind)
     for min = objnorm then (min objnorm min)
     for minind = ind then (if (= objnorm min) ind minind)
     finally (return minind)))
