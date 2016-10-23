package edu.koritsas.slimemold;

import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.shapefile.DirectedIrrigationNetwork;
import edu.koritsas.slimemold.shortestpath.DirecredPhysarumPolycephalumSP;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicGraph;
import org.geotools.graph.traverse.standard.DirectedDijkstraIterator;
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
 * Created by ilias on 13/10/2016.
 */
@RunWith(Parameterized.class)
public class DirectedSlimeSPTest {

    public DirectedSlimeSPTest(DirectedIrrigationNetwork network){
        this.network=network;
    }
    private final Random random=new Random();


    private DirectedIrrigationNetwork network;

    private final  double Io=2;
    private final double γ=1.8;
    private final int numberOfIterations=30000;
    private AbstractDirectedPhysarumPolycephalum slime;
    private Graph graph;
    private DijkstraShortestPathFinder pf;
    private Node source;
    private Node sink;
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
        DirectedIrrigationNetwork network1 = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P1.shp");

        DirectedIrrigationNetwork network2 = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H2.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W2.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P2.shp");

        DirectedIrrigationNetwork network3 = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H3.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W3.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P3.shp");

        DirectedIrrigationNetwork network4 = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P4.shp");
        ////Graph g1 = network1.getGraph();
        // Graph g2 = network2.getGraph();
        // Graph g3 = network3.getGraph();


        return Arrays.asList(new DirectedIrrigationNetwork[]{network1, network2, network3,network4});
    };



    @Before
    public void setUp() throws Exception {

        DirectedGraph graph=network.getBasicGraph();
        List<Node> v=network.getWaterSource();
        source=network.getWaterSource().get(0);;
        sink=network.getHydrants().get(random.nextInt(network.getHydrants().size()));

        slime = new DirecredPhysarumPolycephalumSP(graph,source,sink,Io,γ,numberOfIterations) {
            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry geometry = (Geometry) f.getDefaultGeometry();
                return geometry.getLength();
            }
        };


        DirectedDijkstraIterator.EdgeWeighter weighter = new DirectedDijkstraIterator.EdgeWeighter() {
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