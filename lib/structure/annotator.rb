module Structure
  class Annotator
    def initialize(*fields, &block)
      @fields = fields.flatten
      @block = &block
    end

    def annotate(isoform, residue)
      @block.call isoform, residue
    end
  end
end
