(asdf:defsystem :cl-nsga2
  :author "Ben Diedrich"
  :license "MIT"
  :description "Non Dominated Sorting Genetic Algorithm II (NSGA-II) in Common Lisp"
  :serial t
  :depends-on ()
  :components
  ((:file "package")
   (:file "nsga2")
   (:file "constraint")
   (:file "test")
   (:file "plot")))
