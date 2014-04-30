module Structure
  class Annotator
    attr_accessor :fields, :organism
    def initialize(*fields, &block)
      @fields = fields.flatten
      @block = block
      class << self; self end.send(:define_method, :annotate, &block)
    end
  end
end
