import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.SlimeMold;
import edu.koritsas.slimemold.SlimeMoldSP;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.graph.build.GraphBuilder;
import org.geotools.graph.build.GraphGenerator;
import org.geotools.graph.build.basic.BasicGraphBuilder;
import org.geotools.graph.build.basic.BasicGraphGenerator;
import org.geotools.graph.build.line.BasicLineGraphGenerator;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Graphable;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicGraph;
import org.geotools.graph.traverse.standard.DijkstraIterator;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.List;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

        IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
        //IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P1.shp");
        Graph graph=null;
        try {
            graph =network.getGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }

        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();
        //sinkNodes=sinkNodes.subList(7,8);

        Node source =sourceNodes.get(0);
        Node sink =sinkNodes.get(7);

        SlimeMoldSP slimeMold = new SlimeMoldSP(graph,source,sink,2,1.8,500);
       // SlimeMold slimeMold = new SlimeMold(graph,sourceNodes.get(0),sinkNodes,2,1.8,5000);
        slimeMold.execute();

        slimeMold.showFlowDiagram();
        slimeMold.showConductivityMap();

        try {
            GraphUtils.visualizeGraph(graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(slimeMold.getSolutionCost());


        DijkstraIterator.EdgeWeighter weighter = new DijkstraIterator.EdgeWeighter() {
            @Override
            public double getWeight(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry geom = (Geometry) f.getDefaultGeometry();
                return geom.getLength();
            }
        };

        DijkstraShortestPathFinder pf = new DijkstraShortestPathFinder(graph,source,weighter);
        pf.calculate();

        Path path=pf.getPath(sink);

        BasicGraph g =new BasicGraph(null,path.getEdges());

        try {
            GraphUtils.visualizeGraph(g);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(pf.getCost(sink));

    }
}
