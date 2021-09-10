from shortestPaths import *
import time

def main():
    # APQ Test Code

    # --------------------------
    # pq = AdaptablePQ()
    # pq.add(14, "bed")
    # pq.add(22, "dog")
    # pq.add(18, "fox")
    # pq.add(35, "ant")
    # pq.add(27, "egg")
    # pq.add(24, "cat")
    # # pq.updateCost("ant", 13)
    # print(pq.min().value)
    # # print(pq._heap.left(pq._heap._array[0]).value, pq._heap.right(pq._heap._array[0]).value)
    # print(pq.remove_min().value)
    # print(pq.remove_min().value)
    # print(pq.remove_min().value)
    # print(pq.remove_min().value)
    # print(pq.remove_min().value)
    # --------------------------

    # ENDS

    # Dijkstra test code

    # --------------------------
    # graph = Graph()
    # x=graph.add_vertex("x")
    # y=graph.add_vertex("y")
    # z=graph.add_vertex("z")
    # a=graph.add_vertex("a")
    # b=graph.add_vertex("b")
    # graph.add_edge(x, y, "xy")
    # graph.add_edge(y, z, "yz")
    # graph.add_edge(z, x, "zx")
    # graph.add_edge(z, a, "za")
    # graph.add_edge(y, b, "yb")
    # print(graph.runDijkstra(z))

    # g = graphReaderCorkCity("corkCityData.txt")
    # g.getShortestPath(1669466540, 1147697924, printResults=True, label=True)
    # --------------------------

    # ENDS


    # Cork city route map test code

    # --------------------------

    testCorkCity(timeData=True)

def testCorkCity(timeData=False):
    ids = {}

    ids['wgb'] = 1669466540
    ids['turnerscross'] = 348809726
    ids['neptune'] = 1147697924
    ids['cuh'] = 860206013
    ids['oldoak'] = 358357
    ids['gaol'] = 3777201945
    ids['mahonpoint'] = 330068634

    if timeData:
        startTime = time.gmtime()
        startSecs = time.time()
    # Part 2 method - uses static method, so must reference RouteMap
    graph = RouteMap.graphReaderCorkCity("corkCityData.txt")
    # Part 1 method - much slower due to linear search for getting vertex reference
    # graph=Graph.graphReaderFormatOne("corkCityData.txt")
    graph.sp(ids["mahonpoint"], ids['neptune'], label=True, printPath=True)
    if timeData:
        endTime = time.gmtime()
        endSecs = time.time()
        print("**** Time data ****\nStarted: " + time.asctime(startTime) + "\nEnded: " + time.asctime(endTime) + "\nTotal time (secs): " + str(endSecs-startSecs))




if __name__ == "__main__":
    main()
