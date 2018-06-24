from django.conf.urls import url, include
from rest_framework import routers
from . import rest_views

router = routers.DefaultRouter()
router.register('drugs', rest_views.DrugViewSet)
router.register('drugbankid', rest_views.DrugBankIDViewSet)
router.register('pathway', rest_views.PathwayViewSet)
router.register('target', rest_views.TargetViewSet)
router.register('target-related-mol-chembl', rest_views.ChEMBL_small_molecule_all_infoViewSet)
router.register('target-related-mol-structure', rest_views.ChEMBL_small_moleculeViewSet)

urlpatterns = router.urls
urlpatterns += [
    # url(r'^', include(router.urls)),
    url(r'^target-prediction/', rest_views.target_pred),
    url(r'^bulk-target-prediction', rest_views.bulk_target_pred_2),
    url(r'^feedback/', rest_views.FeedbackCreateView.as_view())
]
