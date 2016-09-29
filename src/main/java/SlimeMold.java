import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.PhysarumPolycephalum;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

import java.util.List;

/**
 * Created by ilias on 29/9/2016.
 */
public class SlimeMold extends PhysarumPolycephalum {

    public SlimeMold(Graph graph, List<Node> sourceNodes, List<Node> sinkNodes, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNodes, sinkNodes, Io, γ, numberOfIterations);
    }

    @Override
    public double getEdgeWeight(Edge e) {
        SimpleFeature f = (SimpleFeature) e.getObject();
        Geometry geometry = (Geometry) f.getDefaultGeometry();
        return geometry.getLength();
    }

    @Override
    public double getEdgeConstraint(Edge e) {
        return 0;
    }
}
