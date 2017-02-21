require "benchmark"
require "kdtree"
require "tempfile"
require "minitest/autorun"
require 'byebug'

#
# create a tree
#

class KdtreeTest < Minitest::Test
  #TMP = "#{Dir.tmpdir}/kdtree_test"

  # def test_node
  #   node = KDTreeNode.new(1.2, 2.3, 15)
  #   assert( (node.x - 1.2).abs < 0.00001 )
  #   assert( (node.y - 2.3).abs < 0.00001 )
  #   assert(node.id == 15)
  # end

  def setup
    @points = (0...1000).map { |i| KDTreeNode.new(rand_coord, rand_coord, i) }
    @kdtree = KDTree.new(@points)
  end

  def find_nearest(points, pt, cnt=1)
    sorted = points.sort_by { |i| distance(i, pt) }
    if cnt == 1
      sorted.first
    else
      sorted.take(cnt)
    end
  end

  def test_nearest
    100.times do
      pt = OpenStruct.new(x: rand_coord, y: rand_coord)

      # One result
      kdpt = @kdtree.nearest(pt.x, pt.y)
      sortpt = find_nearest(@points, pt)
      assert(kdpt.id == sortpt.id, "kdtree nearest did not return the closest result")

      # N results
      kdpts = @kdtree.nearest(pt.x, pt.y, 10)
      sortpts = find_nearest(@points, pt, 10)
      assert(kdpts.map(&:id).sort == sortpts.map(&:id).sort, "kdtree nearest n did not return the closest results")


      # With limit
      kdpts = @kdtree.nearest(pt.x, pt.y, 10, 0.5)
      limited = @points.select{ |p| distance(p, pt) < 0.5 }
      sortpts = find_nearest(limited, pt, 10)
      assert(kdpts.map(&:id).sort == sortpts.map(&:id).sort, "kdtree nearest n with limit did not return the closest results")
    end
  end

  def test_add
    added = []
    1000.times do |i|
      pt = KDTreeNode.new(rand_coord, rand_coord, @points.size+i)
      added << pt
      @kdtree << pt
    end
    100.times do
      pt = OpenStruct.new(x: rand_coord, y: rand_coord)

      # kdtree search
      kdpt = @kdtree.nearest(pt.x, pt.y)

      # slow search
      sortpt = (@points + added).sort_by { |i| distance(i, pt) }.first
      # assert
      kdd = distance(kdpt, pt)
      sortd = distance(sortpt, pt)
      assert((kdd - sortd).abs < 0.0000001, "kdtree didn't return the closest result")
    end
  end

  def test_remove
    removed = []
    500.times do |i|
      node = @points.sample
      @points = @points.delete_if{ |n| n.id == node.id }
      @kdtree = @kdtree.delete(node)
      removed << node
    end
    100.times do
      pt = OpenStruct.new(x: rand_coord, y: rand_coord)

      # kdtree search
      kdpt = @kdtree.nearest(pt.x, pt.y)

      # slow search
      sortpt = @points.sort_by { |i| distance(i, pt) }.first
      # assert
      kdd = distance(kdpt, pt)
      sortd = distance(sortpt, pt)
      assert((kdd - sortd).abs < 0.0000001, "kdtree didn't return the closest result")
    end
  end

  def dont_test_speed
    sizes = [1, 100, 1000, 10000, 100000, 1000000]
    ks = [1, 5, 50, 255]
    sizes.each do |s|
      points = (0...s).map { |i| KDTreeNode.new(rand_coord, rand_coord, i) }

      # build
      Benchmark.bm(17) do |bm|
        kdtree = nil
        bm.report "build" do
          kdtree = KDTree.new(points)
        end

        ks.each do |k|
          bm.report "100 queries (#{k})" do
            total = count = 0
            100.times do
              tm = Time.now
              kdtree.nearest(rand_coord, rand_coord, k, -1)
            end
          end
        end
      end
      puts
    end
  end

  protected

  def distance(a, b)
    x, y = a.x - b.x, a.y - b.y
    x * x + y * y
  end

  def rand_coord
    rand(0) * 10 - 5
  end
end

# running dont_test_speed on my i5 2.8ghz:
#
#                         user     system      total        real
# build               3.350000   0.020000   3.370000 (  3.520528)
# persist             0.150000   0.020000   0.170000 (  0.301963)
# read                0.280000   0.000000   0.280000 (  0.432676)
# 100 queries (1)     0.000000   0.000000   0.000000 (  0.000319)
# 100 queries (5)     0.000000   0.000000   0.000000 (  0.000412)
# 100 queries (50)    0.000000   0.000000   0.000000 (  0.001417)
# 100 queries (255)   0.000000   0.000000   0.000000 (  0.006268)
