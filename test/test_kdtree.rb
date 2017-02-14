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

  # def setup
  #   @points = [
  #     KDTreeNode.new(6, 3, 0),
  #     KDTreeNode.new(2, 5, 1),
  #     KDTreeNode.new(3, 8, 2),
  #     KDTreeNode.new(8, 9, 3),
  #   ]
  #   @kdtree = KDTree.new(@points)
  # end

  # def test_nearest
  #   pt = OpenStruct.new(x: 7, y: 4)
  #   kdpts = @kdtree.nearest(pt.x, pt.y, 10, -1)
  #   kdpt = kdpts.first
  #   sortpt = @points.sort_by { |i| distance(i, pt) }.first
  #   byebug
  #   kdd = distance(kdpt, pt)
  #   sortd = distance(sortpt, pt)
  #   sortpt
  # end

  # def teardown
  #   File.unlink(TMP) if File.exists?(TMP)
  # end

  def test_nearest
    100.times do
      pt = OpenStruct.new(x: rand_coord, y: rand_coord)

      # kdtree search
      kdpts = @kdtree.nearest(pt.x, pt.y, 10, -1)
      kdpt = kdpts.first

      # slow search
      sortpt = @points.sort_by { |i| distance(i, pt) }.first

      # assert
      kdd = distance(kdpt, pt)
      sortd = distance(sortpt, pt)
      assert((kdd - sortd).abs < 0.0000001, "kdtree didn't return the closest result")
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
      kdpts = @kdtree.nearest(pt.x, pt.y, 10, -1)
      kdpt = kdpts.first

      # slow search
      sortpt = (@points + added).sort_by { |i| distance(i, pt) }.first
      # assert
      kdd = distance(kdpt, pt)
      sortd = distance(sortpt, pt)
      assert((kdd - sortd).abs < 0.0000001, "kdtree didn't return the closest result")
    end
  end

  # def test_tmp
  #   @points = [
  #     KDTreeNode.new(6, 3, 0),
  #     KDTreeNode.new(2, 5, 1),
  #     KDTreeNode.new(3, 8, 2),
  #     KDTreeNode.new(8, 9, 3),
  #     KDTreeNode.new(1, 10, 4),
  #   ]
  #   @kdtree = KDTree.new(@points)

  #   byebug
  #   tree
  # end

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
      kdpts = @kdtree.nearest(pt.x, pt.y, 10, -1)
      kdpt = kdpts.first

      # slow search
      sortpt = @points.sort_by { |i| distance(i, pt) }.first
      # assert
      kdd = distance(kdpt, pt)
      sortd = distance(sortpt, pt)
      assert((kdd - sortd).abs < 0.0000001, "kdtree didn't return the closest result")
    end
  end

  # def test_nearestk
  #   100.times do
  #     pt = [rand_coord, rand_coord]

  #     # kdtree search
  #     list = @kdtree.nearestk(pt[0], pt[1], 5)
  #     kdpt = @points[list.last]

  #     # slow search
  #     sortpt = @points.sort_by { |i| distance(i, pt) }[list.length - 1]

  #     # assert
  #     kdd = distance(kdpt, pt)
  #     sortd = distance(sortpt, pt)
  #     assert((kdd - sortd).abs < 0.0000001, "kdtree didn't return the closest result")
  #   end
  # end

  # def test_persist
  #   # write
  #   File.open(TMP, "w") { |f| @kdtree.persist(f) }
  #   # read
  #   kdtree2 = File.open(TMP, "r") { |f| Kdtree.new(f) }

  #   # now test some random points
  #   100.times do
  #     pt = [rand_coord, rand_coord]
  #     id1 = @kdtree.nearest(*pt)
  #     id2 = kdtree2.nearest(*pt)
  #     assert(id1 == id2, "kdtree2 differed from kdtree")
  #   end
  # end

  # def test_bad_magic
  #   File.open(TMP, "w") { |f| f.puts "That ain't right" }
  #   assert_raises RuntimeError do
  #     File.open(TMP, "r") { |f| Kdtree.new(f) }
  #   end
  # end

  # def test_eof
  #   File.open(TMP, "w") { |f| @kdtree.persist(f) }
  #   bytes = File.read(TMP)

  #   [2, 10, 100].each do |len|
  #     File.open(TMP, "w") { |f| f.write(bytes[0, len]) }
  #     assert_raises EOFError do
  #       File.open(TMP, "r") { |f| Kdtree.new(f) }
  #     end
  #   end
  # end

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