# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{bio-emboss_six_frame_nucleotide_sequences}
  s.version = "0.1.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Ben J Woodcroft"]
  s.date = %q{2011-04-03}
  s.description = %q{a method to get the nucleotide sequence of translations done by the EMBOSS bioinformatics package program transeq.}
  s.email = %q{gmail.com after donttrustben}
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.rdoc"
  ]
  s.files = [
    ".document",
    "Gemfile",
    "LICENSE.txt",
    "README.rdoc",
    "Rakefile",
    "VERSION",
    "lib/bio-emboss_six_frame_nucleotide_sequences.rb",
    "lib/bio/sequence/emboss_six_frame_nucleotide_sequences.rb",
    "test/data/test.fa",
    "test/data/test_transeq_6frame.fa",
    "test/helper.rb",
    "test/test_bio-emboss_six_frame_nucleotide_sequences.rb"
  ]
  s.homepage = %q{http://github.com/wwood/bioruby-emboss_six_frame_nucleotide_sequences}
  s.licenses = ["MIT"]
  s.require_paths = ["lib"]
  s.rubygems_version = %q{1.6.2}
  s.summary = %q{a method to get the nucleotide sequence of translations done by the EMBOSS bioinformatics package program transeq.}
  s.test_files = [
    "test/helper.rb",
    "test/test_bio-emboss_six_frame_nucleotide_sequences.rb"
  ]

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<bio>, [">= 1.4.1"])
      s.add_development_dependency(%q<shoulda>, [">= 0"])
      s.add_development_dependency(%q<bundler>, ["~> 1.0.0"])
      s.add_development_dependency(%q<jeweler>, ["~> 1.5.2"])
      s.add_development_dependency(%q<rcov>, [">= 0"])
    else
      s.add_dependency(%q<bio>, [">= 1.4.1"])
      s.add_dependency(%q<shoulda>, [">= 0"])
      s.add_dependency(%q<bundler>, ["~> 1.0.0"])
      s.add_dependency(%q<jeweler>, ["~> 1.5.2"])
      s.add_dependency(%q<rcov>, [">= 0"])
    end
  else
    s.add_dependency(%q<bio>, [">= 1.4.1"])
    s.add_dependency(%q<shoulda>, [">= 0"])
    s.add_dependency(%q<bundler>, ["~> 1.0.0"])
    s.add_dependency(%q<jeweler>, ["~> 1.5.2"])
    s.add_dependency(%q<rcov>, [">= 0"])
  end
end

