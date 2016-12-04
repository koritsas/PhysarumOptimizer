package edu.koritsas.slimemold;

import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumSP;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicGraph;
import org.geotools.graph.traverse.standard.DijkstraIterator;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Random;

/**
 * Created by ilias on 11/10/2016.
 */
@RunWith(Parameterized.class)
public class SlimeMoldSPTest {
    public SlimeMoldSPTest(IrrigationNetwork network) throws IOException {
       this.network=network;


    }
    private final Random random=new Random();


    private IrrigationNetwork network;

    private final  double Io=1;
    private final double γ=1.8;
    private final int numberOfIterations=300;
    private PhysarumPolycephalumSP slime;
    private Graph graph;
    private DijkstraShortestPathFinder pf;
    private Node source;
    private Node sink;
    private double absoluteThreshold =1E-10;
    private double relativeThreshold =1E-10;
    @After
    public void tearDown() throws Exception {
        network=null;
        pf=null;
        slime=null;
        source=null;
        sink=null;

    }

    @Parameterized.Parameters
    public static Collection graphs() throws IOException {
        IrrigationNetwork network1 = new IrrigationNetwork("ParametrizedTests/H1.shp", "ParametrizedTests/W1.shp", "ParametrizedTests/P1.shp");

        IrrigationNetwork network2 = new IrrigationNetwork("ParametrizedTests/H2.shp", "ParametrizedTests/W2.shp", "ParametrizedTests/P2.shp");

        IrrigationNetwork network3 = new IrrigationNetwork("ParametrizedTests/H3.shp", "ParametrizedTests/W3.shp", "ParametrizedTests/P3.shp");

        IrrigationNetwork network4 = new IrrigationNetwork("ParametrizedTests/H4.shp", "ParametrizedTests/W4.shp", "ParametrizedTests/P4.shp");
        ////Graph g1 = network1.getGraph();
       // Graph g2 = network2.getGraph();
       // Graph g3 = network3.getGraph();


        return Arrays.asList(new IrrigationNetwork[]{network1, network2, network3,network4});
    };



    @Before
    public void setUp() throws Exception {

        Graph graph=network.getBasicGraph();
        List<Node> v=network.getWaterSource();
       source=network.getWaterSource().get(0);;
       sink=network.getHydrants().get(random.nextInt(network.getHydrants().size()));

        slime = new PhysarumPolycephalumSP(graph,source,sink,absoluteThreshold,relativeThreshold,numberOfIterations) {
            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry geometry = (Geometry) f.getDefaultGeometry();
                return geometry.getLength();
            }
        };


        DijkstraIterator.EdgeWeighter weighter = new DijkstraIterator.EdgeWeighter() {
            @Override
            public double getWeight(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry geometry = (Geometry) f.getDefaultGeometry();
                return geometry.getLength();
            }
        };

        pf = new DijkstraShortestPathFinder(graph,source,weighter);

    }



    @Test
    public void execute() throws Exception {





        pf.calculate();

        Path path=pf.getPath(sink);




        System.out.println(pf.getCost(sink));

        slime.execute();
         Graph slimeGraph=slime.getGraph();
        System.out.println(slime.getSolutionCost());
        BasicGraph dg = new BasicGraph(null,path.getEdges());

        Collection<Edge> pathEdges =path.getEdges();

        Assert.assertTrue(pathEdges.containsAll(slimeGraph.getEdges()));
        Assert.assertTrue(slimeGraph.getEdges().containsAll(pathEdges));

    }

}