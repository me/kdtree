Gem::Specification.new do |s|
  s.name        = "kdtree"
  s.version     = "0.1"

  s.authors     = ["Ivan Pirlik"]
  s.email       = ["ivan.pirlik@deliveroo.co.uk"]
  s.homepage    = "http://github.com/me/kdtree"
  s.license     = "MIT"
  s.summary     = "2d kdtree with insert and delete."
  s.description = <<EOF
EOF

  s.rubyforge_project = "kdtree"
  s.add_development_dependency "minitest"
  s.add_development_dependency "rake-compiler"

  s.files      = `git ls-files`.split("\n")
  s.test_files = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.extensions = ["ext/kdtree/extconf.rb"]
  s.require_paths = ["lib"]
end
