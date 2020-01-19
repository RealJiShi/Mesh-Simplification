# Mesh Simplification

MeshLab source code flow

- `src\meshlabplugins\filter_meshing\meshfilter.cpp:995`

  - 参数：
    - `OptimalPlacement`
    - `NormalCheck`: `NormalThrRad = PI/4.0;`
    - `QualityThr`

- `QuadricSimplification(m.cm, TargetFaceNum, lastq_Selected, pp, cb);`

  - `m.cm`: mesh
  - `TargetFaceNum`: target number
  - `pp`: 参数

  ```c++
  // Create LocalOptimization
  vcg::LocalOptimization<CMeshO> DeciSession(m, &pp);
  // Init
  DeciSession.Init<tri::MyTriEdgeCollapse>();
  // Optimization
  while (DeciSession.DoOptimization() && m.fn > TargetFaceNum)
  {
      cb(100 - 100 * (m.fn - TargetFaceNum) / (faceToDel), "Simplifying...");
  };
  // Finalize
  DeciSession.Finalize<tri::MyTriEdgeCollapse>();
  ```

- `DeciSession.Init<tri::MyTriEdgeCollapse>();`

  - 实际调用`complex\algorithms\local_optimization\tri_edge_collapse_quadric.h:226`
  - 注意这里的数据结构，存储应为一个half-edge模型，其中
    - `if ((x.V0() < x.V1())`目的是为了判断winding，只添加v0发射出的向量

  ```c++
  static void Init(TriMeshType &m, HeapType &h_ret, BaseParameterClass *_pp)
  {
      // Find the list of faces sharing the vertex
      vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
      // Get face border for the vertex face
      vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);
      // quadric matrix
      InitQuadric(m, pp);
      // create and add heap element
      h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(), x.V1()), TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::GlobalMark(), _pp)));
  }
  ```

  其中的`HeapElem`实质上是一个`TriEdgeCollapseQuadric`(包含一个pair以及当前pair的各种运算)：

  ```c++
  inline TriEdgeCollapseQuadric(const VertexPair &p, int i, BaseParameterClass *pp)
  {
      this->localMark = i;
      this->pos = p;
      this->_priority = ComputePriority(pp);
  }
  ```

  先看看每个`HeapElem`的优先级是如何计算的(注意质量与法向正确性系数)：

  ```c++
  ScalarType ComputePriority(BaseParameterClass *_pp)
  {
      // first store all faces' normal for a pair(v0, v1)
      // store pair ...
      
      // compute position for this pair
      CoordType newPos = ComputePosition(_pp);
      v[0]->P() = v[1]->P() = newPos;
      
      // for all faces in v0(except v1)
      // compute the dot between the new face and the old face, 
      // compute the quality of the new face(quality = minEdgeLength/triangleArea);
      // ...
      
      // for all faces in v1(except v0)
      // compute the dot between the new face and the old face, 
      // compute the quality of the new face(quality = minEdgeLength/triangleArea);
      // ...
      
      // get minQuality and minCos
      // ...
      
      // get the quadric error with function: vQv
      QuadErr = std::max(QuadErr, pp->QuadricEpsilon);
      if (QuadErr <= pp->QuadricEpsilon)
      {
          QuadErr = -1 / Distance(OldPos0, OldPos1);
      }
      
      // finally, priority can be calculated:
      error = (ScalarType)(QuadErr / (MinQual * MinCos));
      
      // restore the v[0] and v[1]
      this->_priority = error;
      return this->_priority;
  }
  ```

  - 再看看`HeapElem`

    ```c++
    struct HeapElem
    {
        inline HeapElem() { locModPtr = NULL; }
        ~HeapElem() {}
        // pointer to instance of local modifier
        // as TriEdgeCollapseQuadric
        LocModPtrType locModPtr;
        float pri;
        inline HeapElem(LocModPtrType _locModPtr)
        {
            locModPtr = _locModPtr;
            pri = float(locModPtr->Priority());
        };
        /// STL heap has the largest element as the first one.
        /// usually we mean priority as an error so we should invert the comparison
        inline bool operator<(const HeapElem &h) const
        {
            return (pri > h.pri);
            //return (locModPtr->Priority() < h.locModPtr->Priority());
        }
        bool IsUpToDate() const
        {
            return locModPtr->IsUpToDate();
        }
    };
    ```

- `DeciSession.DoOptimization()`

  - `complex\algorithms\local_optimization.h:245`

  ```c++
  bool DoOptimization()
  {
      start = clock();
      nPerfmormedOps = 0;
      while (!GoalReached() && !h.empty())
      {
          if (h.size() > m.SimplexNumber() * HeapSimplexRatio)
              ClearHeap();
          std::pop_heap(h.begin(), h.end());
          LocModPtrType locMod = h.back().locModPtr;
          currMetric = h.back().pri;
          h.pop_back();
  
          if (locMod->IsUpToDate())
          {
              // check if it is feasible
              if (locMod->IsFeasible(this->pp))
              {
                  nPerfmormedOps++;
                  locMod->Execute(m, this->pp);
                  locMod->UpdateHeap(h, this->pp);
              }
           }
           //else printf("popped out unfeasible\n");
           delete locMod;
      }
      return !(h.empty());
  }
  ```

  - 首先，如果heap的数量比要简化的数量还多，那么清空heap，随之simplification终止
  - `locMod`在这里就是一个`TriEdgeCollapseQuadric`，表示一个pair
  - 取heap中priority最高的`TriEdgeCollapseQuadric`，即Quadric Error最小的pair进行更新。

- `complex\algorithms\local_optimization\tri_edge_collapse_quadric.h:191`

  ```c++
  void Execute(TriMeshType &m, BaseParameterClass *_pp)
  {
      QH::Qd(this->pos.V(1)) += QH::Qd(this->pos.V(0));
      // v0 is deleted and v1 take the new position
      EdgeCollapser<TriMeshType, VertexPair>::Do(m, this->pos, ComputePosition(_pp)); 
  }
  ```

  - 首先，更新Quadric Matrix, $Q_{1} = Q_{0} + Q_{1}$
  - 随后计算最优位置，并进行`edge collapse`

- `complex\algorithms\local_optimization\tri_edge_collapse_quadric.h:177`

  ```c++
  CoordType ComputePosition(BaseParameterClass *_pp)
  {
      QParameter *pp = (QParameter *)_pp;
      CoordType newPos = (this->pos.V(0)->P() + this->pos.V(1)->P()) / 2.0;
      if ((QH::Qd(this->pos.V(0)).Apply(newPos) + QH::Qd(this->pos.V(1)).Apply(newPos)) > 200.0 * pp->QuadricEpsilon)
          newPos = ComputeMinimal();
      return newPos;
  }
  ```

  先取`v[0]`与`v[1]`的中点，计算`newPos`的误差，如果误差稍大，则可以解方程：
  $$
  V = Q^{-1}*[0\ 0\ 0\ 1]^{T}
  $$
  代码如下：

  ```c++
  CoordType ComputeMinimal()
  {
      typename TriMeshType::VertexType *v[2];
      v[0] = this->pos.V(0);
      v[1] = this->pos.V(1);
      QuadricType q = QH::Qd(v[0]);
      q += QH::Qd(v[1]);
      Point3<QuadricType::ScalarType> x;
      bool rt = q.Minimum(x);
      if (!rt)
      { // if the computation of the minimum fails we choose between the two edge points and the middle one.
          Point3<QuadricType::ScalarType> x0 = Point3d::Construct(v[0]->P());
          Point3<QuadricType::ScalarType> x1 = Point3d::Construct(v[1]->P());
          x.Import((v[0]->P() + v[1]->P()) / 2);
          double qvx = q.Apply(x);
          double qv0 = q.Apply(x0);
          double qv1 = q.Apply(x1);
          if (qv0 < qvx)
              x = x0;
          if (qv1 < qvx && qv1 < qv0)
              x = x1;
      }
      return CoordType::Construct(x);
  }
  ```

  - 计算最优解
  - 如果最优解计算失败，比较`v0`、`v1`和`(v0+v1)/2`的Quadric error，找出其中的最小值最为新的点。

- `complex\algorithms\edge_collapse.h:258`

  ```c++
  static int Do(TriMeshType &m, VertexPair &c, const Point3<ScalarType> &p)
  {
      EdgeSet es;
      FindSets(c, es);
      typename VFIVec::iterator i;
      int n_face_del = 0;
      for (i = es.AV01().begin(); i != es.AV01().end(); ++i)
      {
          FaceType &f = *((*i).f);
          assert(f.V((*i).z) == c.V(0));
          vcg::face::VFDetach(f, ((*i).z + 1) % 3);
          vcg::face::VFDetach(f, ((*i).z + 2) % 3);
          Allocator<TriMeshType>::DeleteFace(m, f);
          n_face_del++;
      }
      //set Vertex Face topology
      for (i = es.AV0().begin(); i != es.AV0().end(); ++i)
      {
          (*i).f->V((*i).z) = c.V(1);
          (*i).f->VFp((*i).z) = (*i).f->V((*i).z)->VFp();
          (*i).f->VFi((*i).z) = (*i).f->V((*i).z)->VFi();
          (*i).f->V((*i).z)->VFp() = (*i).f;
          (*i).f->V((*i).z)->VFi() = (*i).z;
      }
      Allocator<TriMeshType>::DeleteVertex(m, *(c.V(0)));
      c.V(1)->P() = p;
      return n_face_del;
  }
  ```

  - `FindSets(c, es);`找出对应的edge：
    - `VFIVec &AV0()`: Face in v0(except v1).
    - `VFIVec &AV1()`: Face in v1(except v0).
    - `VFIVec &AV01()`: Face between v0 and v1.
  - 拿到所有edge之后，将`AV01`的三角面全部删除
  - 最终，将`AV0`中的三角面的v0索引改成v1
  - 删除顶点v0，更新v1的位置信息（newPos）

- `complex\algorithms\local_optimization\tri_edge_collapse_quadric.h:451`

  ```c++
  inline void UpdateHeap(HeapType &h_ret, BaseParameterClass *_pp)
  {
      this->GlobalMark()++;
      VertexType *v[2];
      v[0] = this->pos.V(0);
      v[1] = this->pos.V(1);
      v[1]->IMark() = this->GlobalMark();
      // First loop around the surviving vertex to unmark the Visit flags
      for (VFIterator vfi(v[1]); !vfi.End(); ++vfi)
      {
          vfi.V1()->ClearV();
          vfi.V2()->ClearV();
      }
      // Second Loop
      for (VFIterator vfi(v[1]); !vfi.End(); ++vfi)
      {
          if (!(vfi.V1()->IsV()) && vfi.V1()->IsRW())
          {
              vfi.V1()->SetV();
              AddCollapseToHeap(h_ret, vfi.V0(), vfi.V1(), _pp);
          }
          if (!(vfi.V2()->IsV()) && vfi.V2()->IsRW())
          {
              vfi.V2()->SetV();
              AddCollapseToHeap(h_ret, vfi.V2(), vfi.V0(), _pp);
          }
          if (vfi.V1()->IsRW() && vfi.V2()->IsRW())
              AddCollapseToHeap(h_ret, vfi.V1(), vfi.V2(), _pp);
      } // end second loop around surviving vertex.
  }
  ```

  - 遍历v1的所有的三角形，判断三角形的每一条边是否符合条件，如符合将其加入Heap中：

    ```c++
    inline void AddCollapseToHeap(HeapType &h_ret, VertexType *v0, VertexType *v1, BaseParameterClass *_pp)
    {
        QParameter *pp = (QParameter *)_pp;
        h_ret.push_back(HeapElem(new MYTYPE(VertexPair(v0, v1), this->GlobalMark(), _pp)));
        std::push_heap(h_ret.begin(), h_ret.end());
    }
    ```

    